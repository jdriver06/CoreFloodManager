import numpy as np
from PyQt5.QtWidgets import QMainWindow, QVBoxLayout, QHBoxLayout, QFrame, QPushButton, QWidget, \
    QLineEdit, QDialog, QComboBox, QLabel, QGroupBox, QTabWidget, QListWidget, QGridLayout, QCheckBox, QMessageBox, \
    QTextEdit, QScrollArea
from PyQt5.Qt import QIcon, Qt, QFont
from PyQt5.QtGui import QMouseEvent, QCloseEvent, QKeyEvent
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtChart import QChart, QChartView, QLineSeries, QValueAxis
import fluid
import flood
import core
import reports
import j_utils as utils
from numpy import isnan, array, round, append, cumsum, where, isinf, divide, pi, sqrt
# from docx import Document
# from docx.shared import Inches


__version__ = '0.1.0'


class InjectionFluidsView(QDialog):

    def __init__(self, slm: utils.SignalListManager, parent=None):
        super(InjectionFluidsView, self).__init__(parent)

        self.setWindowTitle(' ')

        cls = [[fluid.FluidView], [fluid.InjectionFluid]]
        args = [[], []]
        self.list_widget = utils.SignalListManagerWidget(self, slm, ['Injection Fluids'], cls, *args)

        pmv = self.parent().parent()
        signal_lists = [pmv.brine_list.signal_lists[0]]
        for s_list in pmv.chemical_list.signal_lists:
            signal_lists.append(s_list)
        signal_lists.append(pmv.oil_list.signal_lists[1])
        signal_lists.append(pmv.oil_list.signal_lists[2])
        signal_lists.append(pmv.gas_list.signal_lists[0])
        self.list_widget.pm_list.signal_lists[0].set_ref_lists(signal_lists)

        layout = QVBoxLayout()
        layout.addWidget(self.list_widget)
        self.setLayout(layout)
        if parent is not None:
            sf = parent.sf
            # sf = parent.width() / 1000.
        else:
            sf = 1.
        self.setFixedSize(int(sf * 280), int(sf * 350))
        self.sf = sf

        self.show()

    def closeEvent(self, _):

        try:
            self.parent().if_view = None
            for icon in self.parent().flood_icons:
                icon.reset_if_combobox()
        except Exception as e:
            print(e)


class MultiFloodsView(QDialog):

    def __init__(self, slm: utils.SignalListManager, parent=None):
        super(MultiFloodsView, self).__init__(parent=parent)
        self.setWindowTitle(' ')

        cls = [[flood.MultiCoreFloodBuilder], [flood.MultiCoreFlood]]
        args = []
        self.list_widget = utils.SignalListManagerWidget(self, slm, ['Multi-Floods'], cls, *args)

        layout = QVBoxLayout()
        layout.addWidget(self.list_widget)
        self.setLayout(layout)

        self.show()

    def closeEvent(self, _):

        try:
            self.parent().mf_view = None
        except Exception as e:
            print(e)


class FloodExperiment:

    def __init__(self, name: str, c: core.Core, pm=None):

        self.name = name
        self.objective = ''
        self.notes = ''
        self.project_manager_list = pm
        self.view_class = FloodExperimentView
        self.core = c
        self.petro_parameters = {'C': [2, 2, None], 'perm': [1000., 1000., None], 'krw0': [0.25, 0.25, None],
                                 'kro0': [0.5, 0.5, None], 'Swr': [0.2, 0.2, None], 'Sor': [0.2, 0.2, None],
                                 'SoI': [0.2, 0.2, None], 'nw': [2.5, 2.5, None], 'no': [2.5, 2.5, None],
                                 'phi': [0.25, 0.25, None], 'IPV': [0.0, 0.0, None]}

        tLphi = 0.
        for section in c.core_sections:
            tLphi += section.length * section.est_porosity

        phi_est = tLphi / c.my_length()
        self.petro_parameters['phi'] = [phi_est, phi_est, None]
        self.initial_state = 0
        self.floods = []
        self.ref_objects = []
        self.floods_list = utils.SignalListManager()
        self.floods_list.add_signal_list()
        self.injection_fluids_list = utils.SignalListManager()
        self.injection_fluids_list.add_signal_list()
        self.multi_floods_list = utils.SignalListManager()
        self.multi_floods_list.add_signal_list()
        self.multi_floods_list.signal_lists[0].set_ref_lists([self.floods_list.signal_lists[0]])
        self.ift_dict = dict()

    def find_multiflood(self, f: flood.CoreFlood):
        """ This function returns the multi-flood claiming flood f if it exists, otherwise, None is returned. """

        for mf in self.multi_floods_list.signal_lists[0].objects:
            if f in mf.ref_floods:
                return mf

        return None

    def binned_floods(self, f: flood.CoreFlood=None) -> list:
        """ This function returns a list of floods and multi-floods where multi-floods replace their constituent floods
        unless the flood f is in that multi-flood. As an example, consider a flood_experiment with floods f1, f2, f3,
        f4, f5, and f6, as well as multi-floods mf1(f1, f2) and mf2(f5, f6). A call to binned_floods(f5) will return a
        list [mf1, f3, f4, f5, f6], while a call binned_floods(f3) would return [mf1, f3, f4, mf2]. """

        # If there are no multi-floods, just return the list of floods.
        if not self.multi_floods_list.signal_lists[0].objects:
            return self.floods

        # If the flood f is not a multi-flood, get a handle to any multi-flood associated with it (may be None even if
        # the "if" statement evaluates as True).
        if f is not None and not isinstance(f, flood.MultiCoreFlood):
            f_mf = self.find_multiflood(f)
        else:
            f_mf = None

        # Initialize a list of binned floods.
        bfs = []

        # Build the binned list of floods. A series of base floods in a multi-flood should be replaced by a single
        # handle to the multi-flood unless the flood f is in that series.
        for fld in self.floods:
            mf = self.find_multiflood(fld)
            # If fld is not part of a multi-flood, it should be appended. It should also be appended if it is part of a
            # multi-flood associated with the flood "f". Else, append the multi-flood mf associated with the flood fld
            # unless it has already been appended.
            if mf is None or mf == f_mf:
                bfs.append(fld)
            else:
                if mf not in bfs:
                    bfs.append(mf)

        return bfs

    def get_most_recent_oil_if(self, f: flood.CoreFlood) -> fluid.InjectionFluid:

        oil_if = None

        bf = self.binned_floods(f)
        i = bf.index(f)

        for fld in reversed(bf[:i + 1]):
            flu = fld.get_last_fluid()
            if isinstance(flu.specific_fluid, fluid.OilInjectionFluid):
                oil_if = flu
                break

        return oil_if

    def get_ift(self, f: flood.CoreFlood, flu: fluid.InjectionFluid = None) -> float:

        if not isinstance(f, flood.MultiCoreFlood):
            flu = f.fluid

        if isinstance(flu.specific_fluid, fluid.OilInjectionFluid):
            return np.nan

        oil_if = self.get_most_recent_oil_if(f)

        if oil_if is None:
            return np.nan

        if (oil_if, flu) in self.ift_dict.keys():
            return self.ift_dict[(oil_if, flu)]

        if isinstance(f, flood.MultiCoreFlood):
            fluids_list = f.get_fluids_list()
            min_ift = np.nan
            for fluidi in fluids_list:
                if fluidi == flu:
                    continue
                ifti = np.nan
                if (oil_if, fluidi) in self.ift_dict.keys():
                    ifti = self.ift_dict[(oil_if, fluidi)]
                if isnan(min_ift) or ifti < min_ift:
                    min_ift = ifti

            if not isnan(min_ift):
                return min_ift

        return oil_if.specific_fluid.oil_sample.get_ift()

    def get_pv(self) -> float:

        return self.petro_parameters['phi'][1] * self.core.my_bulk_volume()

    def calculate_saturation(self, flood_list: list, i: int, j: int, sw0: float, direction: bool) -> float:
        """ The function calculate_saturation calculates the saturation at the ith flood from reference flood j of
        saturation sw0 incrementing in the forward direction (True) or reverse direction (False). """

        sw = sw0
        pv = self.get_pv()

        if direction:
            # If we are incrementing forward from the reference flood j, take a slice of all floods after j and up to
            # and including the flood in question (i)
            floods = flood_list[j + 1:i + 1]
        else:
            # In the reverse direction, take everything after the flood in question (i) up through the reference flood.
            # The saturation of each flood is the saturation AT THE END of the flood. Reverse the elements of the slice
            # so that they may be incremented from the reference flood to the flood in question.
            floods = flood_list[i + 1: j + 1]
            floods = floods[::-1]

        print('floods: ', floods, floods[0].fluid.specific_fluid)

        for fld in floods:
            # print(fld.name, 'Sw start:', sw)
            if isinstance(fld.fluid.specific_fluid, fluid.BrineInjectionFluid):
                if direction:
                    # print('adding saturation', sw)
                    # print('oil produced: ', fld.effluent.oil_produced())
                    sw += fld.effluent.oil_produced() / pv
                else:
                    # print('subtracting saturation', sw)
                    # print('oil produced: ', fld.effluent.oil_produced())
                    sw -= fld.effluent.oil_produced() / pv
            else:
                if direction:
                    # print('subtracting saturation', sw)
                    # print('brine produced: ', fld.effluent.brine_produced())
                    sw -= fld.effluent.brine_produced() / pv
                else:
                    # print('adding saturation', sw)
                    # print('brine produced: ', fld.effluent.brine_produced())
                    sw += fld.effluent.brine_produced() / pv
            # print(fld.name, 'Sw end:', sw)
        print('Sw = {:.3f}'.format(sw))
        return sw

    def get_flood_saturation(self, f: flood.CoreFlood):
        """The function get_flood_saturation returns a flood's aqueous saturation."""

        flood_list = self.binned_floods(f)

        if f not in flood_list:
            return -1.
        else:
            print('flood {} is in list'.format(f.name))

        if f.n_phases == 1:
            # If the flood is single-phase, return a saturation of 1 for a brine injection fluid and 0 for oil.
            if isinstance(f.fluid.specific_fluid, fluid.BrineInjectionFluid):
                print('hey')
                return 1.
            else:
                return 0.

        i = flood_list.index(f)
        # Reference floods (refs) will contain flood index (col 1) and saturation (col 2)
        refs = [[], []]

        if self.initial_state != 1:
            # An initial state of "1" would signify two-phase. Otherwise, it is "0" (aqueous) or "2" (oil).
            refs[0].append(0)
            if self.initial_state == 0:
                refs[1].append(1.)
            else:
                refs[1].append(0.)

        phi = self.petro_parameters['phi'][1]

        for j, fld in enumerate(flood_list):
            # Scan through floods and fill in saturations. A flood's pore_volume attribute refers to pore volume from
            # R.I. or aqueous tracer.
            if not isnan(fld.pore_volume) and fld.use_for_ref:
                # A flood must have a defined pore volume to be a reference.
                refs[0].append(j)
                # print(fld.pore_volume, phi, self.core.my_bulk_volume())
                refs[1].append(fld.pore_volume / (phi * self.core.my_bulk_volume()))

        if i in refs[0]:
            # If the flood we want saturation information for is a reference, then just return that saturation.
            return refs[1][refs[0].index(i)]

        if not refs[1]:
            # If there are no reference floods/saturations, calculate from estimate of SwI.
            swI = 1. - self.petro_parameters['SoI'][1]

            if i == 0:
                return swI
            print('no refs, calculating saturation: ', i, swI, True)
            return self.calculate_saturation(flood_list, i, 0, swI, True)

        if len(refs[1]) == 1:
            # If there's one reference, find the saturation starting from that reference.
            if refs[0][0] < i:
                direction = True

            else:
                direction = False

            print('one ref, calculating saturation:', i, refs[0][0], refs[1][0], direction)
            return self.calculate_saturation(flood_list, i, refs[0][0], refs[1][0], direction)

        else:
            # If there's more than one reference, try to use both.
            r0 = array(refs[0])
            r1 = array(refs[1])

            # Construct lists of reference floods before and after the flood being interrogated.
            greater = r0[r0 > i]
            lesser = r0[r0 < i]
            g = greater.tolist()
            l = lesser.tolist()

            # Initialize the "less than" (forward-propagated) and "greater than" (backward propagated) saturations to a
            # flag saturation of -1.
            swl = -1.
            swg = -1.

            if l:
                # If there are prior floods that are references, select the corresponding slice of saturations and pass
                # the nearest (last) reference index and saturation to the calculate saturation function, which returns
                # a forward-propagated saturation.
                r1l = r1[r0 < i]
                print('prior ref, calculating saturation:', i, l[-1], r1l[-1], True)
                swl = self.calculate_saturation(flood_list, i, l[-1], r1l[-1], True)

            if g:
                # If there are subsequent floods that are references, pass the first flood index and saturation to the
                # calculate_saturation function.
                r1g = r1[r0 > i]
                print('subsequent ref, calculating saturation:', i, g[0], r1g[0], False)
                swg = self.calculate_saturation(flood_list, i, g[0], r1g[0], False)

            if swl > 0.:
                # If there was a prior reference flood, use it with the subsequent reference flood, if there is one.
                # Otherwise, use only the subsequent flood.
                if swg > 0.:
                    sw = (swl + swg) / 2.
                else:
                    sw = swl
            else:
                sw = swg

            # print('sw', sw)
            return sw

    def get_initial_saturation(self) -> float:

        fl0 = self.floods[0]
        sw1 = self.get_flood_saturation(fl0)
        pv = self.get_pv()

        if isinstance(fl0.fluid.specific_fluid, fluid.BrineInjectionFluid):
            return sw1 - fl0.effluent.oil_produced() / pv

        return sw1 + fl0.effluent.brine_produced() / pv

    def get_core_mass(self, cut_num: int) -> float:

        mass = 0.
        phi = self.petro_parameters['phi'][1]
        for cs in self.core.core_sections:
            mass += (1. - phi) * cs.my_bulk_volume(cut_num) * core.CoreSection.densities[cs.lithology]

        return mass

    def load_from_pickle(self, pickle):

        self.core = pickle.core
        self.floods = []

        for fss in pickle.flood_save_structures:
            self.floods.append(flood.CoreFlood(fss.name, fss.fluid, fss.core))
            self.floods[-1].load_from_pickle(fss)


class FloodExperimentSaveStructure:

    def __init__(self, fe: FloodExperiment):
        self.core = core
        self.flood_save_structures = []
        for f in fe.floods:
            self.flood_save_structures.append(flood.CoreFloodSaveStructure(f))


class FloodExperimentView(QMainWindow):

    def __init__(self, fe: FloodExperiment, parent=None):
        super(FloodExperimentView, self).__init__(parent=parent)

        icon = QIcon('Coreholder Cropped.jpg')
        self.setWindowIcon(icon)
        self.setWindowTitle('Flood Experiment: ' + fe.name + ' on core ' + fe.core.name)
        self.setFixedSize(1000, 400)

        sf = 1.
        if parent is not None:
            self.move(parent.x(), parent.y())
            sf = parent.width() / 900.
            self.setFixedSize(int(sf * 1000), int(sf * 400))

        self.sf = sf

        c_widget = QWidget(parent=self)
        c_layout = QVBoxLayout()
        c_widget.setLayout(c_layout)

        f_widget = QWidget(parent=self)
        f_layout = QHBoxLayout()
        f_widget.setLayout(f_layout)

        b_widget = QWidget(parent=self)
        b_widget.setFixedHeight(int(sf * 50))
        b_layout = QHBoxLayout()
        b_layout.setContentsMargins(0, 0, 0, 0)
        b_widget.setLayout(b_layout)

        c_layout.addWidget(f_widget)
        c_layout.addWidget(b_widget)
        self.setCentralWidget(c_widget)

        self.left_button = QPushButton('<')
        self.left_button.setFixedWidth(20)
        self.left_button.clicked.connect(self.left_button_click)
        self.right_button = QPushButton('>')
        self.right_button.setFixedWidth(20)
        self.right_button.clicked.connect(self.right_button_click)
        self.frame = QFrame()
        self.frame.setStyleSheet('background: white')
        self.frame.setLayout(QHBoxLayout())

        f_layout.addWidget(self.left_button)
        f_layout.addWidget(self.frame)
        f_layout.addWidget(self.right_button)

        self.insert_button = QPushButton('Insert')
        self.insert_button.clicked.connect(self.insert_flood_wrapper)
        self.delete_button = QPushButton('Delete')
        self.delete_button.clicked.connect(self.delete_flood)
        self.multi_floods_button = QPushButton('MultiFloods')
        self.multi_floods_button.clicked.connect(self.multi_floods_view)
        self.injection_fluids_button = QPushButton('Injection Fluids')
        self.injection_fluids_button.clicked.connect(self.injection_fluids_view)
        self.petrophysics_button = QPushButton('Petrophysics')
        self.petrophysics_button.clicked.connect(self.petrophysics_view)
        self.plan_button = QPushButton('Plan')
        self.plan_button.clicked.connect(self.plan_view)
        self.report_button = QPushButton('Summary')
        self.report_button.clicked.connect(self.summary_view)
        self.export_button = QPushButton('Export')
        self.export_button.clicked.connect(self.export_to_excel_wrapper)
        font = self.insert_button.font()
        font.setPointSize(10)
        self.insert_button.setFont(font)
        self.delete_button.setFont(font)
        self.multi_floods_button.setFont(font)
        self.injection_fluids_button.setFont(font)
        self.petrophysics_button.setFont(font)
        self.plan_button.setFont(font)
        self.report_button.setFont(font)
        self.export_button.setFont(font)
        b_layout.addWidget(self.insert_button)
        b_layout.addWidget(self.delete_button)
        b_layout.addWidget(self.multi_floods_button)
        b_layout.addWidget(self.injection_fluids_button)
        b_layout.addWidget(self.petrophysics_button)
        b_layout.addWidget(self.plan_button)
        b_layout.addWidget(self.report_button)
        b_layout.addWidget(self.export_button)
        # self.report_button.setEnabled(False)

        self.core = fe.core
        self.experiment = fe
        self.flood_icons = []
        self.flood_forms = []
        self.petro_view = None
        self.summary_viewer = None
        self.selected_icon = -1
        self.start_icon = -1
        self.total_floods_created = 0

        n_floods = len(fe.floods)

        for i in range(n_floods):
            self.insert_flood_icon(fe.floods[i])

        self.selected_icon = -1

        if fe.floods:
            for i in range(n_floods):
                self.selected_icon = i
                self.flood_icons[i].sync_if_combobox_to_flood()

        self.if_view = None
        self.mf_view = None

        self.control_key_pressed = False

        self.show()

    def set_visible_flood_icons(self):

        if not self.flood_icons:
            return

        for i, icon in enumerate(self.flood_icons):
            if i < self.start_icon or i > self.start_icon + 2:
                icon.setVisible(False)
            else:
                icon.setVisible(True)

    def correct_start_icon(self):

        if len(self.experiment.floods) < 3:
            self.start_icon = 0
        elif self.selected_icon > self.start_icon + 2:
            self.start_icon = self.selected_icon - 2
        elif self.start_icon > len(self.experiment.floods) - 3:
            self.start_icon -= 1

    def right_button_click(self):

        if self.start_icon + 3 < len(self.flood_icons):
            self.start_icon += 1

        self.set_visible_flood_icons()

    def left_button_click(self):

        if self.start_icon > 0:
            self.start_icon -= 1

        self.set_visible_flood_icons()

    def insert_flood_wrapper(self):

        i = self.selected_icon

        if i + 1 < len(self.experiment.floods):
            f1 = self.experiment.floods[i]
            f2 = self.experiment.floods[i + 1]
            mf1 = self.experiment.find_multiflood(f1)
            mf2 = self.experiment.find_multiflood(f2)
            print('insert_flood_wrapper', mf1, mf2)

            if mf1 == mf2 and mf1 is not None:
                QMessageBox(parent=self, text='Cannot insert new flood into Multi-Flood ' + mf1.name).exec_()
                return

        try:
            self.insert_flood(None)
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()

    def insert_flood(self, f=None):

        j = self.selected_icon + 1
        self.total_floods_created += 1

        if f is None:
            f = flood.CoreFlood('flood ' + str(self.total_floods_created),
                                fluid.InjectionFluid('', [fluid.Brine('')], [[], 'brine-based']),
                                self.core, self.experiment)

            self.experiment.floods.insert(j, f)
            if len(self.experiment.floods) > 1:
                self.experiment.floods[-1].manual_offsets[:] = self.experiment.floods[-2].manual_offsets[:]
                self.experiment.floods[-1].cut_num = self.experiment.floods[-2].cut_num
            else:
                experiments = self.experiment.project_manager_list.pm_list.signal_lists[0].objects
                i = experiments.index(self.experiment)

                if i > 0:
                    k = i
                    while k > 0:
                        k -= 1
                        prev_experiment = experiments[k]
                        if prev_experiment.core == self.core and prev_experiment.floods:
                            k = -1
                            self.experiment.floods[0].cut_num = prev_experiment.floods[-1].cut_num

        else:
            self.experiment.floods.insert(j, f)

        self.insert_flood_icon(f)
        self.update_flood_phases()

    def insert_flood_icon(self, f=None):

        if self.selected_icon != -1:
            self.flood_icons[self.selected_icon].select(False)

        j = self.selected_icon + 1

        if f is None:
            name = 'flood ' + str(self.total_floods_created)
        else:
            name = f.name

        self.flood_forms.insert(j, None)
        self.experiment.floods[j].name = name

        flood_icon = FloodIcon(self.experiment.floods[-1], self)
        self.flood_icons.append(flood_icon)
        self.frame.layout().addWidget(flood_icon)

        for i, icon in enumerate(self.flood_icons):
            icon.flood = self.experiment.floods[i]
            icon.cff = self.flood_forms[i]
            icon.name_edit.setText(icon.flood.name)
            icon.temp_edit.setText(str(icon.flood.temperature))
            icon.bp_edit.setText(str(icon.flood.back_pressure))
            icon.v_edit.setText(str(icon.flood.planned_volume))
            kk = 0
            for k, injection_fluid in enumerate(self.experiment.injection_fluids_list.signal_lists[0].objects):
                if injection_fluid == icon.flood.fluid:
                    kk = k + 1
            icon.if_combobox.setCurrentIndex(kk)
            icon.cut_combobox.setCurrentIndex(icon.flood.cut_num)
            icon.rev_checkbox.setChecked(icon.flood.reverse_flow)
            icon.select(False)

        self.flood_icons[j].select()
        self.selected_icon = j

        self.correct_start_icon()
        self.set_visible_flood_icons()

    def delete_flood(self):

        f = self.experiment.floods[self.selected_icon]
        mf = self.experiment.find_multiflood(f)
        if mf is not None:
            QMessageBox(parent=self, text='Cannot delete flood in use by Multi-Flood ' + mf.name)
            return

        self.frame.layout().removeWidget(self.flood_icons[self.selected_icon])
        self.experiment.floods.pop(self.selected_icon)
        icon = self.flood_icons.pop(self.selected_icon)
        self.flood_forms.pop(self.selected_icon)
        icon.close()
        if self.selected_icon > len(self.flood_icons) - 1:
            self.selected_icon = self.selected_icon - 1

        if self.flood_icons:
            self.flood_icons[self.selected_icon].select()

        self.correct_start_icon()
        self.set_visible_flood_icons()
        self.update_flood_phases()

    def update_selected_flood(self, fi, right_click: bool=False):

        i = self.flood_icons.index(fi)

        if self.selected_icon != i:
            self.flood_icons[self.selected_icon].select(False)
            self.selected_icon = i
            self.flood_icons[i].select()
        elif right_click:
            self.flood_icons[self.selected_icon].select(False)
            self.selected_icon = -1

    def update_flood_phases(self):

        fe = self.experiment

        for i, fl in enumerate(fe.floods):

            # print(fl.fluid)

            if not isinstance(fl.fluid, fluid.InjectionFluid):
                if i == 0:
                    if fe.initial_state == 0 or fe.initial_state == 2:
                        fl.n_phases = 1
                    else:
                        fl.n_phases = 2
                else:
                    fl.n_phases = fe.floods[i - 1].n_phases
                continue

            if i == 0:

                if fe.initial_state == 0:
                    if isinstance(fl.fluid.specific_fluid, fluid.BrineInjectionFluid):
                        fl.n_phases = 1
                    elif isinstance(fl.fluid.specific_fluid, fluid.OilInjectionFluid):
                        fl.n_phases = 2
                    # else:
                        # print('fluid is not of either type.')
                        # print(type(fl.fluid.specific_fluid))
                elif fe.initial_state == 2:
                    if isinstance(fl.fluid.specific_fluid, fluid.BrineInjectionFluid):
                        fl.n_phases = 2
                    elif isinstance(fl.fluid.specific_fluid, fluid.OilInjectionFluid):
                        fl.n_phases = 1
                    # else:
                    #     print('fluid is not of either type.')
                    #     print(type(fl.fluid.specific_fluid))
                else:
                    fl.n_phases = 2

            else:

                if fe.floods[i - 1].n_phases == 2:
                    fl.n_phases = 2
                elif not isinstance(fe.floods[i - 1].fluid, fluid.InjectionFluid):
                    print('fluid is not an injection fluid.')
                    j = i - 1
                    while j > 0:
                        j -= 1
                        if isinstance(fe.floods[j].fluid, fluid.InjectionFluid):
                            if isinstance(fl.fluid.specific_fluid, type(fe.floods[j].fluid.specific_fluid)):
                                fl.n_phases = fe.floods[j].n_phases
                            else:
                                fl.n_phases = 2
                            break
                        elif j == 0:
                            fl.n_phases = fe.floods[j].n_phases
                elif isinstance(fl.fluid.specific_fluid, type(fe.floods[i - 1].fluid.specific_fluid)):
                    fl.n_phases = 1
                else:
                    fl.n_phases = 2

        for mf in fe.multi_floods_list.signal_lists[0].objects:
            mf.n_phases = mf.ref_floods[0].n_phases

            # print(i, fl.name, fl.n_phases)
            # print(self.experiment.get_flood_saturation(fl))

    def update_flood_forms_saturation_text(self):

        for form in self.flood_forms:
            if form is not None:
                form.update_saturation_text()

    def injection_fluids_view(self):
        if self.if_view is None:
            self.if_view = InjectionFluidsView(self.experiment.injection_fluids_list, self)

    def multi_floods_view(self):
        if self.mf_view is None:
            self.mf_view = MultiFloodsView(self.experiment.multi_floods_list, self)

    def petrophysics_view(self):

        if self.petro_view is None:
            try:
                self.petro_view = PetrophysicsView(self.experiment, self, self.sf)
            except Exception as e:
                QMessageBox(parent=self, text=str(e)).exec_()

    def plan_view(self):

        file_name = utils.file_open_dlg(self.experiment.name + ' Plan', '(PDF) *.pdf', True)

        try:
            reports.create_plan(self.experiment, file_name)
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()

    def report_view(self):

        builder = ReportBuilder(self)
        builder.exec_()

    def summary_view(self):

        if self.summary_viewer is not None:
            return

        self.summary_viewer = SummaryViewer(self)
        self.summary_viewer.show()

    def export_to_excel_wrapper(self):

        try:
            self.export_to_excel()
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()

    def export_to_excel(self):
        file_name = utils.file_open_dlg('', '.xlsx', True)
        if not file_name:
            return
        import xlsxwriter
        # import numpy as np
        with xlsxwriter.Workbook(file_name + '.xlsx') as workbook:
            flds = self.experiment.binned_floods()
            for fl in flds:
                if fl.p_data[0].size > 1:

                    name = fl.name + ' (data)'
                    if len(name) > 31:
                        name = name[-31:]
                    ws = workbook.add_worksheet(name=name)
                    p_data = fl.p_data
                    t_data = fl.time_data
                    v_data = fl.flow_data
                    frs = fl.flow_rate
                    fr_data = t_data.copy()
                    for tp, fr in frs.items():
                        fr_data[t_data >= tp] = fr

                    is_rf_data = False
                    if fl.rf_data[0].size and fl.rf_data[0].size == fl.p_data[0].size:
                        is_rf_data = True

                    if not is_rf_data:
                        line_0 = ['Raw', '', '', '', '', 'Offset']
                        line_1 = ['Flow Rate [mL/min]', 'Vol. [mL]', 'Vol. [PV]',
                                  'Whole [psi]', 'Sec 1 [psi]', 'Sec 2 [psi]', 'Sec 3 [psi]', 'Sec 4 [psi]',
                                  'Whole [psi]', 'Sec 1 [psi]', 'Sec 2 [psi]', 'Sec 3 [psi]', 'Sec 4 [psi]']
                    else:
                        line_0 = ['Raw', '', '', '', '', 'Offset', '', '', '', '', 'R.F.']
                        line_1 = ['Flow Rate [mL/min]', 'Vol. [mL]', 'Vol. [PV]',
                                  'Whole [psi]', 'Sec 1 [psi]', 'Sec 2 [psi]', 'Sec 3 [psi]', 'Sec 4 [psi]',
                                  'Whole [psi]', 'Sec 1 [psi]', 'Sec 2 [psi]', 'Sec 3 [psi]', 'Sec 4 [psi]',
                                  'Whole', 'Sec 1', 'Sec 2', 'Sec 3', 'Sec 4']

                    ws.write_row(0, 3, line_0)
                    ws.write_row(1, 0, line_1)

                    iv = [fl.get_nn(i) for i in range(5)]

                    write_data = [fr_data, v_data, v_data / fl.get_experiment_total_pore_volume(),
                                  p_data[0], p_data[iv[1]], p_data[iv[2]], p_data[iv[3]], p_data[iv[4]]]

                    for i in range(5):
                        write_data.append(p_data[i] - fl.used_offsets[i])
                    if is_rf_data:
                        for i in range(5):
                            d = fl.rf_data[i].copy()
                            infs = where(isinf(d))
                            d[infs] = 0.
                            write_data.append(d)

                    for col_num, data in enumerate(write_data):
                        ws.write_column(2, col_num, data)

                    series_names = ['Whole', 'Sec1', 'Sec2', 'Sec3', 'Sec4']
                    colors = ['#4F81BD', '#C0504D', '#9BBB59', '#8064A2', '#F79646']

                    name2 = fl.name + ' (Plot)'
                    if len(name2) > 31:
                        name2 = name2[-31:]

                    cht_sht = workbook.add_chartsheet(name=name2)
                    cht = workbook.add_chart({'type': 'scatter'})
                    columns = ['I', 'J', 'K', 'L', 'M']
                    for series_name, c, color in zip(series_names, columns, colors):
                        cht.add_series({'name': series_name,
                                        'categories': '=\'' + name + '\'!$C$3:$C$' + str(v_data.size + 2),
                                        'values': '=\'' + name + '\'!$' + c + '$3:$' + c + '$' + str(v_data.size + 2),
                                        'marker': {'type': 'circle',
                                                   'border': {'color': color},
                                                   'fill': {'color': color}}})

                    cht.set_x_axis({'name': 'Pore Volumes Injected',
                                    'name_font': {'name': 'Arial', 'size': 16},
                                    'num_format': '0.00',
                                    'num_font': {'name': 'Courier New', 'size': 14},
                                    'major_gridlines': {'visible': True, 'line': {'width': 0.75, 'color': '#D9D9D9'}},
                                    'minor_gridlines': {'visible': True, 'line': {'width': 0.75, 'color': '#F2F2F2'}}})
                    cht.set_y_axis({'name': chr(916) + u'P [psi]',
                                    'name_font': {'name': 'Arial', 'size': 16},
                                    'num_format': '0.0',
                                    'num_font': {'name': 'Courier New', 'size': 14},
                                    'major_gridlines': {'visible': True, 'line': {'width': 0.75, 'color': '#D9D9D9'}},
                                    'minor_gridlines': {'visible': True, 'line': {'width': 0.75, 'color': '#F2F2F2'}}})
                    cht.set_legend({'position': 'overlay_right',
                                    'font': {'name': 'Arial', 'size': 12},
                                    'fill': {'color': '#DCE6F2'},
                                    'border': {'none': True}})

                    cht.set_chartarea({'border': {'none': True}})

                    cht_sht.set_chart(cht)

                    if is_rf_data:

                        name2 = fl.name + ' (R.F. Plot)'
                        if len(name2) > 31:
                            name2 = name2[-31:]

                        cht_sht = workbook.add_chartsheet(name=name2)
                        cht = workbook.add_chart({'type': 'scatter'})
                        columns = ['N', 'O', 'P', 'Q', 'R']
                        for series_name, c, color in zip(series_names, columns, colors):
                            cht.add_series({'name': series_name,
                                            'categories': '=\'' + name + '\'!$C$3:$C$' + str(v_data.size + 2),
                                            'values': '=\'' + name + '\'!$' + c + '$3:$' + c + '$' + str(
                                                v_data.size + 2),
                                            'marker': {'type': 'circle',
                                                       'border': {'color': color},
                                                       'fill': {'color': color}}})

                        cht.set_x_axis({'name': 'Pore Volumes Injected',
                                        'name_font': {'name': 'Arial', 'size': 16},
                                        'num_format': '0.00',
                                        'num_font': {'name': 'Courier New', 'size': 14},
                                        'major_gridlines': {'visible': True,
                                                            'line': {'width': 0.75, 'color': '#D9D9D9'}},
                                        'minor_gridlines': {'visible': True,
                                                            'line': {'width': 0.75, 'color': '#F2F2F2'}}})
                        cht.set_y_axis({'name': chr(916) + u'P [psi]',
                                        'name_font': {'name': 'Arial', 'size': 16},
                                        'num_format': '0.0',
                                        'num_font': {'name': 'Courier New', 'size': 14},
                                        'major_gridlines': {'visible': True,
                                                            'line': {'width': 0.75, 'color': '#D9D9D9'}},
                                        'minor_gridlines': {'visible': True,
                                                            'line': {'width': 0.75, 'color': '#F2F2F2'}}})
                        cht.set_legend({'position': 'overlay_right',
                                        'font': {'name': 'Arial', 'size': 12},
                                        'fill': {'color': '#DCE6F2'},
                                        'border': {'none': True}})

                        cht.set_chartarea({'border': {'none': True}})

                        cht_sht.set_chart(cht)

                name = fl.name + ' analysis'
                if len(name) > 31:
                    name = name[-31:]
                ws2 = workbook.add_worksheet(name=name)
                ws2.write_row(0, 0, ['Pressure Offsets'])
                ws2.write_row(1, 0, ['Whole [psi]', 'Sec 1 [psi]', 'Sec 2 [psi]', 'Sec 3 [psi]', 'Sec 4 [psi]'])
                ws2.write_row(2, 0, fl.used_offsets)

                try:
                    sw = self.experiment.get_flood_saturation(fl)
                    ws2.write_row(4, 0, ['Sw at end of flood'])
                    ws2.write_row(5, 0, [sw])
                except Exception as e:
                    print(e)

                if fl.plateau_x[0]:
                    ws2.write_row(7, 0, ['Plateau Data'])

                    if not fl.plateau_whole_rf:
                        ws2.write_row(8, 0, ['Begin [min]', 'End [min]', 'Rate [mL/min]', 'Visc. [cP]',
                                             'Whole [psi]', 'Sec 1 [psi]', 'Sec 2 [psi]', 'Sec 3 [psi]', 'Sec 4 [psi]'])
                    else:
                        ws2.write_row(8, 0, ['Begin [min]', 'End [min]', 'Rate [mL/min]', 'Visc. [cP]',
                                             'Whole [psi]', 'Sec 1 [psi]', 'Sec 2 [psi]', 'Sec 3 [psi]', 'Sec 4 [psi]',
                                             'Whole R.F.', 'Sec 1 R.F.', 'Sec 2 R.F.', 'Sec 3 R.F.', 'Sec 4 R.F.'])

                    rate_at_end = []
                    tps = list(fl.flow_rate.keys())
                    for px in fl.plateau_x[1]:
                        tp0 = tps[0]
                        for tp in tps:
                            if px > tp:
                                tp0 = tp
                        rate_at_end.append(fl.flow_rate[tp0])

                    shear = fl.get_plateau_shear_rates(self.experiment, [])
                    try:
                        fluids = fl.get_plateau_fluids([])
                        fluids[0].get_viscosity(fl.temperature, shear[0])
                    except Exception as e:
                        QMessageBox(parent=self, text=str(e)).exec_()
                        return

                    mus = []
                    for s, f in zip(shear, fluids):
                        try:
                            mus.append(f.get_viscosity(fl.temperature, s))
                        except Exception as e:
                            print(e)
                            mus.append(0.)

                    write_data = [fl.plateau_x[0], fl.plateau_x[1], rate_at_end, mus,
                                  fl.plateau_whole - fl.used_offsets[0], fl.plateau_sec1 - fl.used_offsets[1],
                                  fl.plateau_sec2 - fl.used_offsets[2], fl.plateau_sec3 - fl.used_offsets[3],
                                  fl.plateau_sec4 - fl.used_offsets[4]]
                    if fl.plateau_whole_rf:
                        write_data.append(fl.plateau_whole_rf)
                        write_data.append(fl.plateau_sec1_rf)
                        write_data.append(fl.plateau_sec2_rf)
                        write_data.append(fl.plateau_sec3_rf)
                        write_data.append(fl.plateau_sec4_rf)

                    for col_num, data in enumerate(write_data):
                        ws2.write_column(9, col_num, data)

                if list(fl.effluent.volumes[0]):
                    name = fl.name + ' prod. fluids'
                    if len(name) > 31:
                        name = name[-31:]
                    ws3 = workbook.add_worksheet(name=name)
                    ws3_headers = []
                    ws3_headers.append('PV Cent.')
                    ws3_headers.append('Aq. [mL]')
                    ws3_headers.append('M.I. [mL]')
                    ws3_headers.append('Tot. [mL]')
                    ws3_headers.append('Oil [mL]')
                    ws3_headers.append('Fin. Aq. [mL]')
                    ws3_headers.append('Fin. Tot. [mL]')
                    ws3_headers.append('Fin. Oil [mL]')
                    ws3_headers.append('O.C.')
                    ws3_headers.append('So')
                    ws3_headers.append('C.O.R.')
                    ws3_headers.append('')
                    ws3_headers.append('pH')
                    ws3_headers.append('R.I.')
                    ws3_headers.append('ORP [mV]')
                    ws3_headers.append('[cP]')
                    ws3_headers.append('A520')
                    ws3_headers.append('')
                    if list(fl.effluent.ppm):
                        if fl.effluent.cp_or_abs_used == 0:
                            ws3_headers.append('cp coeffs')
                        else:
                            ws3_headers.append('abs coeffs')
                        ws3_headers.append('Polymer [ppm]')
                        ws3_headers.append('Retention [mcg/g]')
                        ws3_headers.append('')
                    ws3_headers.append('Samp [g]')
                    ws3_headers.append('Dil [g]')
                    ws3_headers.append('Tot [g]')
                    ws3_headers.append('Na+ [ppm]')
                    ws3_headers.append('K+ [ppm]')
                    ws3_headers.append('Mg++ [ppm]')
                    ws3_headers.append('Ca++ [ppm]')
                    ws3_headers.append('Cl- [ppm]')
                    ws3_headers.append('Br- [ppm]')
                    ws3_headers.append('SO4-- [ppm]')
                    ws3.write_row(0, 0, ws3_headers)

                    ws3.write_column(1, 0, fl.effluent.cumulative_volume_centered() /
                                     fl.get_experiment_total_pore_volume())
                    col_count = 1
                    for col, vol_data in enumerate(fl.effluent.volumes):
                        if col == 1:
                            for row, dp in enumerate(vol_data):
                                if not isnan(dp):
                                    rv = fl.effluent.revisions[col, row]
                                    ws3.write_row(row + 1, col + 1, [dp + rv])
                        else:
                            if vol_data.shape[0] == fl.effluent.revisions.shape[1]:
                                net_data = vol_data + fl.effluent.revisions[col, :]
                            else:
                                net_data = vol_data + fl.effluent.revisions[col, :-1]
                            ws3.write_column(1, col + 1, list(net_data))
                        col_count += 1
                    fow = fl.effluent.final_oil_and_water()
                    ws3.write_column(1, col_count + 0, fow[0, :])
                    ws3.write_column(1, col_count + 1, fow[1, :])
                    ws3.write_column(1, col_count + 2, fow[2, :])
                    oc = fl.effluent.oil_cut()
                    for row, dp in enumerate(oc):
                        if not isnan(dp):
                            ws3.write_row(row + 1, col_count + 3, [dp])
                    sw = self.experiment.get_flood_saturation(fl)
                    oil_start = (1. - sw) * fl.get_experiment_total_pore_volume() + fl.effluent.oil_produced()
                    so = (oil_start - fl.effluent.cumulative_oil_volume()) / fl.get_experiment_total_pore_volume()
                    ws3.write_column(1, col_count + 4, so)
                    if oil_start:
                        ws3.write_column(1, col_count + 5, fl.effluent.cumulative_oil_volume() / oil_start)
                    col_count += 7

                    col_count_start = col_count
                    for col, m_data in enumerate(fl.effluent.measurements):
                        for row, dp in enumerate(list(m_data)):
                            if not isnan(dp):
                                ws3.write_row(row + 1, col + col_count_start, [dp])
                        col_count += 1
                    col_count += 1
                    if fl.effluent.cp_or_abs_used == 0 and list(fl.effluent.ppm_to_cp_fit_params):
                        ws3.write_column(1, col_count, fl.effluent.ppm_to_cp_fit_params.tolist())
                        col_count += 1
                    elif list(fl.effluent.ppm_to_abs_fit_params):
                        ws3.write_column(1, col_count, fl.effluent.ppm_to_abs_fit_params.tolist())
                        col_count += 1
                    if list(fl.effluent.ppm):
                        ws3.write_column(1, col_count, fl.effluent.ppm)
                        col_count += 1
                    if not isnan(fl.effluent.retention):
                        ws3.write_column(1, col_count, [fl.effluent.retention])
                        col_count += 2
                    for col, i_data in enumerate(fl.effluent.ions):
                        for row, dp in enumerate(list(i_data)):
                            if not isnan(dp):
                                ws3.write_row(row + 1, col_count, [dp])
                        col_count += 1

                    name2 = fl.name + ' (OR Plot)'
                    if len(name2) > 31:
                        name2 = name2[-31:]

                    cht_sht = workbook.add_chartsheet(name=name2)
                    cht = workbook.add_chart({'type': 'scatter'})
                    columns = ['I', 'J', 'K']
                    series_names = ['OC', 'So', 'COR']
                    colors = ['red', '#9BBB59', '#4F81BD']

                    for series_name, c, color in zip(series_names, columns, colors):
                        cht.add_series({'name': series_name,
                                        'categories': '=\'' + name + '\'!$A$2:$A$' + str(so.size),
                                        'values': '=\'' + name + '\'!$' + c + '$2:$' + c + '$' + str(so.size),
                                        'line': {'color': color},
                                        'marker': {'type': 'none'}})

                    cht.set_x_axis({'name': 'Pore Volumes Injected',
                                    'name_font': {'name': 'Arial', 'size': 16},
                                    'num_format': '0.00',
                                    'num_font': {'name': 'Courier New', 'size': 14},
                                    'major_gridlines': {'visible': True, 'line': {'width': 0.75, 'color': '#D9D9D9'}},
                                    'minor_gridlines': {'visible': True, 'line': {'width': 0.75, 'color': '#F2F2F2'}}})
                    cht.set_y_axis({'name': 'Oil Cut, So, Cumulative Oil Recovery [%]',
                                    'name_font': {'name': 'Arial', 'size': 16},
                                    'num_format': '0.0%',
                                    'num_font': {'name': 'Courier New', 'size': 14},
                                    'major_gridlines': {'visible': True, 'line': {'width': 0.75, 'color': '#D9D9D9'}},
                                    'minor_gridlines': {'visible': True, 'line': {'width': 0.75, 'color': '#F2F2F2'}}})
                    cht.set_legend({'position': 'overlay_right',
                                    'font': {'name': 'Arial', 'size': 12},
                                    'fill': {'color': '#DCE6F2'},
                                    'border': {'none': True}})

                    cht.set_chartarea({'border': {'none': True}})

                    cht_sht.set_chart(cht)

    def keyPressEvent(self, a0: QKeyEvent):

        if a0.key() == Qt.Key_Control:
            self.control_key_pressed = True

    def keyReleaseEvent(self, a0: QKeyEvent):

        if a0.key() == Qt.Key_Control:
            self.control_key_pressed = False

    def closeEvent(self, _):
        for icon in self.flood_icons:
            if icon.cff is not None:
                icon.cff.close()
        if self.mf_view is not None:
            self.mf_view.close()
        if self.if_view is not None:
            self.if_view.close()


class FloodIcon(QFrame):

    def __init__(self, f: flood.CoreFlood, fe_view: FloodExperimentView):
        super(FloodIcon, self).__init__()

        self.flood = f
        self.injection_fluid = None
        self.fe_view = fe_view
        self.cff = None
        self.background_color = 'black'

        sf = fe_view.width() / 1000
        self.setFixedSize(int(sf * 275), int(sf * 275))
        self.setStyleSheet('background: black; border: 2px solid gray')
        self.setLayout(QVBoxLayout())

        self.name_edit = QLineEdit(parent=self)
        self.name_edit.setText(f.name)
        self.name_edit.setAlignment(Qt.AlignHCenter)
        self.name_edit.setStyleSheet('background: white; border: 0 px solid gray')
        self.name_edit.textChanged.connect(self.change_flood_name)
        font = self.name_edit.font()
        font.setPointSize(12)
        self.name_edit.setFont(font)
        self.layout().addWidget(self.name_edit)

        font.setPointSize(11)

        self.if_combobox = QComboBox(parent=self)
        self.if_combobox.setStyleSheet('background: white; border: 0 px solid gray')
        self.reset_if_combobox()
        self.if_combobox.currentIndexChanged.connect(self.injection_fluid_select)
        self.if_combobox.setFont(font)
        self.layout().addWidget(self.if_combobox)

        temp_group = QGroupBox(parent=self)
        temp_group.setStyleSheet('border: 0 px solid gray')
        temp_layout = QHBoxLayout()
        temp_layout.setContentsMargins(0, 0, 0, 0)
        temp_group.setLayout(temp_layout)
        temp_text = QLabel(parent=self, text='T [' + chr(176) + 'C]: ')
        temp_text.setFont(font)
        temp_text.setStyleSheet('background: white; border: 0 px solid gray')
        self.temp_edit = QLineEdit(parent=self)
        self.temp_edit.setText(str(f.temperature))
        self.temp_edit.editingFinished.connect(self.update_flood_temp)
        self.temp_edit.setFont(font)
        self.temp_edit.setStyleSheet('background: white; border: 0 px solid gray')
        temp_layout.addWidget(temp_text)
        temp_layout.addWidget(self.temp_edit)
        self.layout().addWidget(temp_group)

        bp_group = QGroupBox(parent=self)
        bp_group.setStyleSheet('border: 0 px solid gray')
        bp_layout = QHBoxLayout()
        bp_layout.setContentsMargins(0, 0, 0, 0)
        bp_group.setLayout(bp_layout)
        bp_text = QLabel(parent=self, text='BP [psig]: ')
        bp_text.setFont(font)
        bp_text.setStyleSheet('background: white; border: 0 px solid gray')
        self.bp_edit = QLineEdit(parent=self)
        self.bp_edit.setText(str(f.back_pressure))
        self.bp_edit.editingFinished.connect(self.update_flood_bp)
        self.bp_edit.setFont(font)
        self.bp_edit.setStyleSheet('background: white; border: 0 px solid gray')
        bp_layout.addWidget(bp_text)
        bp_layout.addWidget(self.bp_edit)
        self.layout().addWidget(bp_group)

        v_group = QGroupBox(parent=self)
        v_group.setStyleSheet('border: 0 px solid gray')
        v_layout = QHBoxLayout()
        v_layout.setContentsMargins(0, 0, 0, 0)
        v_group.setLayout(v_layout)
        v_text = QLabel(parent=self, text='Vol. (est.) [mL]: ')
        v_text.setFont(font)
        v_text.setStyleSheet('background: white; border: 0 px solid gray')
        self.v_edit = QLineEdit(parent=self)
        self.v_edit.setText(str(f.planned_volume))
        self.v_edit.editingFinished.connect(self.update_flood_v)
        self.v_edit.setFont(font)
        self.v_edit.setStyleSheet('background: white; border: 0 px solid gray')
        v_layout.addWidget(v_text)
        v_layout.addWidget(self.v_edit)
        self.layout().addWidget(v_group)

        cut_group = QGroupBox(parent=self)
        cut_group.setStyleSheet('border: 0 px solid gray')
        cut_layout = QHBoxLayout()
        cut_layout.setContentsMargins(0, 0, 0, 0)
        cut_group.setLayout(cut_layout)
        cut_text = QLabel(parent=self, text='Cut #:')
        cut_text.setFont(font)
        cut_text.setStyleSheet('background: white; border: 0 px solid gray')
        self.cut_combobox = QComboBox(parent=self)
        for i in range(len(f.core.cut_refs)):
            self.cut_combobox.addItem(str(i + 1))
        self.cut_combobox.setFont(font)
        self.cut_combobox.setStyleSheet('background: white; border: 0 px solid gray')
        self.cut_combobox.currentIndexChanged.connect(self.cut_select)
        cut_layout.addWidget(cut_text)
        cut_layout.addWidget(self.cut_combobox)
        self.layout().addWidget(cut_group)

        # rev_group = QGroupBox(parent=self)
        # rev_layout = QHBoxLayout()
        # rev_group.setLayout(rev_layout)
        # rev_text = QLabel(parent=self, text='Reverse Flow?')
        # rev_text.setFont(font)
        # rev_text.setStyleSheet('background: white')
        self.rev_checkbox = QCheckBox(parent=self, text='Reverse Flow?')
        self.rev_checkbox.setChecked(self.flood.reverse_flow)
        self.rev_checkbox.setFont(font)
        self.rev_checkbox.setStyleSheet('background: white; border: 0 px solid gray')
        self.rev_checkbox.clicked.connect(self.rev_flow)
        # rev_layout.addWidget(rev_text)
        # rev_layout.addWidget(self.rev_checkbox)
        self.layout().addWidget(self.rev_checkbox)

    def reset_if_combobox(self):
        if_list = list(' ')
        i_names = self.fe_view.experiment.injection_fluids_list.signal_lists[0].item_names
        for i_n in i_names:
            if_list.append(i_n)
        self.if_combobox.clear()
        self.if_combobox.addItems(if_list)

        found = False
        if self.injection_fluid is not None:
            for i, i_n in enumerate(i_names):
                if self.injection_fluid.name == i_n:
                    self.if_combobox.setCurrentIndex(i+1)
                    found = True

        if not found:
            self.injection_fluid = None

    def sync_if_combobox_to_flood(self):

        objs = self.fe_view.experiment.injection_fluids_list.signal_lists[0].objects

        for i, obj in enumerate(objs):
            if obj == self.flood.fluid:
                self.if_combobox.setCurrentIndex(i+1)

    def change_flood_name(self):
        self.flood.name = self.name_edit.text()

    def update_flood_temp(self):

        prev_temp = self.flood.temperature
        temp = self.temp_edit.text()

        try:
            temp = float(temp)
        except ValueError:
            self.temp_edit.setText(str(prev_temp))
            return

        if self.flood.effluent.read_temperature == self.flood.temperature:
            self.flood.effluent.read_temperature = temp
        self.flood.temperature = temp

        mf = self.fe_view.experiment.find_multiflood(self.flood)
        if mf is not None:
            if self.flood == mf.ref_floods[0]:
                mf.temperature = temp

    def update_flood_bp(self):

        prev_bp = self.flood.back_pressure
        bp = self.bp_edit.text()

        try:
            bp = float(bp)
        except ValueError:
            self.bp_edit.setText(str(prev_bp))
            return

        self.flood.back_pressure = bp

        mf = self.fe_view.experiment.find_multiflood(self.flood)
        if mf is not None:
            if self.flood == mf.ref_floods[0]:
                mf.back_pressure = bp

        if self.injection_fluid.get_fluid_type() == fluid.FluidType.LIVE_OIL:
            oil = self.injection_fluid.specific_fluid.oil_sample.ref_objects[0]
            if bp < oil.bubble_point:
                self.bp_edit.setStyleSheet('background: yellow')
            else:
                self.bp_edit.setStyleSheet('background: white')

    def update_flood_v(self):

        prev_v = self.flood.planned_volume
        v = self.v_edit.text()

        try:
            v = float(v)
        except ValueError:
            self.v_edit.setText(str(prev_v))
            return

        self.flood.planned_volume = v

    def select(self, on_off: bool=True):
        if on_off:
            self.setStyleSheet('background: ' + self.background_color + '; border: 3px solid yellow')
        else:
            self.setStyleSheet('background: ' + self.background_color + '; border: 3px solid gray')

    def mousePressEvent(self, a0: QMouseEvent):
        super(FloodIcon, self).mousePressEvent(a0)

        if a0.button() == 1:
            bfs = self.fe_view.experiment.binned_floods(self.flood)
            bf_names = []
            for bf in bfs:
                bf_names.append(bf.name)
            print(bf_names)
            self.fe_view.update_selected_flood(self, False)
        if a0.button() == 2:
            self.fe_view.update_selected_flood(self, True)

    def mouseDoubleClickEvent(self, a0: QMouseEvent):
        super(FloodIcon, self).mouseDoubleClickEvent(a0)

        fe_view = self.fe_view
        if fe_view.flood_forms[fe_view.selected_icon] is not None:
            return

        if self.if_combobox.currentIndex() == 0:
            return

        c_flood = fe_view.experiment.floods[fe_view.selected_icon]
        mcf = fe_view.experiment.find_multiflood(c_flood)

        if mcf is None or fe_view.control_key_pressed:
            c_flood.fluid = self.injection_fluid
            try:
                cff = flood.CoreFloodForm(self.injection_fluid, fe_view.core, fe_view, fe_view, c_flood, self)
                fe_view.flood_forms[fe_view.selected_icon] = cff
                cff.show()
                self.cff = cff
            except Exception as e:
                QMessageBox(parent=self, text=str(e)).exec_()

            if fe_view.control_key_pressed:
                fe_view.control_key_pressed = False

        else:

            for fld, ff in zip(fe_view.experiment.floods, fe_view.flood_forms):
                if ff is not None and fld in mcf.ref_floods:
                    return

            print('building cff.')
            try:
                cff = flood.CoreFloodForm(self.injection_fluid, fe_view.core, fe_view, fe_view, mcf, self)
            except Exception as e:
                QMessageBox(parent=self.parent(), text=str(e)).exec_()
                return
            print('cff built.')

            for i, fld in enumerate(fe_view.experiment.floods):
                if fld in mcf.ref_floods:
                    fe_view.flood_forms[i] = cff
                    fe_view.flood_icons[i].cff = cff
            cff.show()

    def injection_fluid_select(self):

        i = self.if_combobox.currentIndex()
        objs = self.fe_view.experiment.injection_fluids_list.signal_lists[0].objects

        if i > 0:
            if objs:
                self.injection_fluid = objs[i-1]
                self.flood.fluid = self.injection_fluid

                if isinstance(self.injection_fluid.specific_fluid, fluid.BrineInjectionFluid):
                    self.background_color = 'lightblue'
                    for additive in self.injection_fluid.ref_objects:
                        if isinstance(additive, fluid.Polymer):
                            self.background_color = 'lightgray'
                        elif isinstance(additive, fluid.Surfactant) or isinstance(additive, fluid.Formulation):
                            self.background_color = '#CDCD6B'
                            break

                elif isinstance(self.injection_fluid.specific_fluid, fluid.OilInjectionFluid):
                    self.background_color = '#783C00'
                    for additive in self.injection_fluid.ref_objects:
                        if isinstance(additive, fluid.Gas):
                            self.background_color = '#B45A00'
                            oil = self.injection_fluid.specific_fluid.oil_sample.ref_objects[0]
                            print(oil.bubble_point, self.flood.back_pressure)
                            if oil.bubble_point > self.flood.back_pressure:
                                self.bp_edit.setStyleSheet('background: yellow')
                            else:
                                self.bp_edit.setStyleSheet('background: white')
                            break
                        else:
                            self.bp_edit.setStyleSheet('background: white')
        else:
            self.background_color = 'black'
            if self.cff is not None:
                print('closing core flood form')
                self.cff.close()

        if self.fe_view.flood_icons[self.fe_view.selected_icon] == self:
            self.setStyleSheet('background: ' + self.background_color + '; border: 3px solid yellow')
        else:
            self.setStyleSheet('background: ' + self.background_color + '; border: 3px solid gray')

        self.fe_view.update_flood_phases()

    def cut_select(self):

        self.flood.cut_num = self.cut_combobox.currentIndex()

    def rev_flow(self):

        self.flood.set_reverse_flow(self.rev_checkbox.isChecked())

        j = self.fe_view.experiment.floods.index(self.flood)
        if self.fe_view.flood_forms[j] is not None:
            self.fe_view.flood_forms[j].update_from_source()
        mf = self.fe_view.experiment.find_multiflood(self.flood)
        if mf is not None:
            for i, icon in enumerate(self.fe_view.flood_icons):
                fl = self.fe_view.experiment.floods[i]
                if fl in mf.ref_floods:
                    fl.reverse_flow = self.rev_checkbox.isChecked()
                    icon.rev_checkbox.setChecked(self.rev_checkbox.isChecked())
            mf.reverse_flow = self.flood.reverse_flow

    def closeEvent(self, a0: QCloseEvent):
        print('closing flood icon', self.cff)

        if self.cff is not None:
            print('closing core flood form')
            self.cff.flood_icon = None
            self.cff.closeEvent(a0)


class ReportBuilder(QDialog):

    def __init__(self, parent: FloodExperimentView):
        super(ReportBuilder, self).__init__(parent=parent)
        self.setFixedSize(200, 200)
        self.setWindowTitle('Report Builder')
        self.submit_report = QPushButton(parent=self, text='Submit')
        lyt = QGridLayout()
        self.setLayout(lyt)
        lyt.addWidget(self.submit_report)
        # self.submit_report.clicked.connect(self.submit_report_clicked)

        self.rec_plot = RecoveryPlot(self, parent.experiment)

    # def submit_report_clicked(self):
    #
    #     QApplication.screens()[0].grabWindow(self.rec_plot.winId()).save('oil_rec.jpg', 'jpg')
    #     document = Document()
    #
    #     document.add_heading(self.parent().experiment.name + ' Report', 0)
    #
    #     p = document.add_paragraph('A plain paragraph having some ')
    #     p.add_run('bold').bold = True
    #     p.add_run(' and some ')
    #     p.add_run('italic.').italic = True
    #
    #     document.add_heading('Heading, level 1', level=1)
    #     document.add_paragraph('Intense quote', style='Intense Quote')
    #
    #     document.add_paragraph(
    #         'first item in unordered list', style='List Bullet'
    #     )
    #     document.add_paragraph(
    #         'first item in ordered list', style='List Number'
    #     )
    #
    #     document.add_picture('oil_rec.jpg', width=Inches(4.5))
    #
    #     records = (
    #         (3, '101', 'Spam'),
    #         (7, '422', 'Eggs'),
    #         (4, '631', 'Spam, spam, eggs, and spam')
    #     )
    #
    #     table = document.add_table(rows=1, cols=3)
    #     hdr_cells = table.rows[0].cells
    #     hdr_cells[0].text = 'Qty'
    #     hdr_cells[1].text = 'Id'
    #     hdr_cells[2].text = 'Desc'
    #     for qty, id, desc in records:
    #         row_cells = table.add_row().cells
    #         row_cells[0].text = str(qty)
    #         row_cells[1].text = id
    #         row_cells[2].text = desc
    #
    #     document.add_page_break()
    #
    #     document.save(self.parent().experiment.name + ' Report.docx')


class RecoveryPlot(QMainWindow):

    def __init__(self, parent: ReportBuilder, experiment: FloodExperiment):
        super(RecoveryPlot, self).__init__(parent=parent)

        self.experiment = experiment

        self.setFixedSize(800, 500)
        self.setWindowTitle('Oil Recovery Chart for ' + self.experiment.name)
        self.chart = QChart(flags=Qt.WindowFlags())
        self.view = QChartView(self.chart)
        self.setCentralWidget(self.view)

        axis_title_font = QFont()
        axis_title_font.setFamily('Courier')
        axis_title_font.setPointSize(14)
        axis_title_font.setBold(True)
        axis_font = QFont()
        axis_font.setFamily('Courier')
        axis_font.setPointSize(11)
        self.x_axis = QValueAxis()
        self.x_axis.setTitleText('Pore Volumes Injected')
        self.x_axis.setTitleFont(axis_title_font)
        self.x_axis.setLabelsFont(axis_font)
        self.y_axis = QValueAxis()
        self.y_axis.setTitleText('%')
        self.y_axis.setTitleFont(axis_title_font)
        self.y_axis.setLabelsFont(axis_font)
        self.chart.setAxisX(self.x_axis)
        self.chart.setAxisY(self.y_axis)
        self.chart.legend().hide()

        expt = self.experiment
        pv = expt.core.my_area() * expt.core.my_length() * expt.petro_parameters['phi'][1]
        cv = 0.
        oil_lines = []

        soi = 0.
        for fld in expt.floods:
            if isinstance(fld.fluid.specific_fluid, fluid.OilInjectionFluid):
                soi = 1. - expt.get_flood_saturation(fld)

        so = array([])
        rec = array([])
        vols = array([])
        last_so = soi * pv
        so_line = QLineSeries()
        rec_line = QLineSeries()

        for fld in expt.floods:

            eff = fld.effluent
            if eff.oil_produced() > 0.:
                oil_line = QLineSeries()
                oil_lines.append(oil_line)
                cum_vols = eff.cumulative_volume()
                cum_oil = cumsum(eff.volumes[2, :])
                so = append(so, last_so - cum_oil)

                vols = append(vols, (cum_vols + cv) / pv)
                last_so = so[-1]
                oil_line.append(utils.series_to_polyline((cum_vols + cv) / pv, eff.oil_cut()))
                # so_line.append(utils.series_to_polyline((cum_vols + cv) / pv, so))
                if cv > 0.:
                    oil_lines[-2].append(utils.series_to_polyline(([cum_vols[0]] + cv) / pv, [eff.oil_cut()[0]]))
                cv += cum_vols[-1]
                self.chart.addSeries(oil_line)
                oil_line.attachAxis(self.x_axis)
                oil_line.attachAxis(self.y_axis)
                pen = oil_line.pen()
                pen.setColor(Qt.red)
                oil_line.setPen(pen)

        rec = append(rec, 1. - divide(so, soi * pv))
        so_line.append(utils.series_to_polyline(vols, so / pv))
        self.chart.addSeries(so_line)
        so_line.attachAxis(self.x_axis)
        so_line.attachAxis(self.y_axis)
        pen = so_line.pen()
        pen.setColor(Qt.darkGreen)
        so_line.setPen(pen)
        rec_line.append(utils.series_to_polyline(vols, rec))
        self.chart.addSeries(rec_line)
        rec_line.attachAxis(self.x_axis)
        rec_line.attachAxis(self.y_axis)
        pen.setColor(Qt.blue)
        rec_line.setPen(pen)
        self.x_axis.setMax(cv / pv)

        self.show()


class PorosityTab(QWidget):

    def __init__(self, parent):
        super(PorosityTab, self).__init__(parent=parent)


class FloodListTab(QWidget):

    list_item_selected = pyqtSignal()

    def __init__(self, parent):
        super(FloodListTab, self).__init__(parent=parent)

        layout = QGridLayout()
        self.setLayout(layout)

        layout.addWidget(QLabel(parent=self, text='Floods:'), 0, 0)
        self.flood_list = QListWidget(parent=self)
        self.flood_list.currentItemChanged.connect(self.list_item_select)
        layout.addWidget(self.flood_list, 1, 0, 1, 2)

    def list_item_select(self):

        self.list_item_selected.emit()


class PermeabilityTab(FloodListTab):

    def __init__(self, parent):
        super(PermeabilityTab, self).__init__(parent=parent)

        self.list_item_selected.connect(self.display_list_item_data)
        self.offset_combobox = QComboBox(parent=self)
        self.offset_combobox.setEnabled(False)
        self.permeability_combobox = QComboBox(parent=self)
        self.permeability_combobox.setEnabled(False)
        self.rf_combobox = QComboBox(parent=self)
        self.rf_combobox.setEnabled(False)
        self.permeability_text = QLabel(parent=self, text='')
        self.r_permeability_text = QLabel(parent=self, text='')
        fl_list = parent.parent().experiment.binned_floods()
        for fl in fl_list:
            self.flood_list.addItem(fl.name)
            self.offset_combobox.addItem(fl.name)
            self.permeability_combobox.addItem(fl.name)
            self.rf_combobox.addItem(fl.name)

        self.layout().addWidget(QLabel(parent=self, text='Pressure Offset Source: '), 2, 0)
        self.layout().addWidget(self.offset_combobox, 2, 1)
        self.layout().addWidget(QLabel(parent=self, text='Permeability Source: '), 3, 0)
        self.layout().addWidget(self.permeability_combobox, 3, 1)
        self.layout().addWidget(QLabel(parent=self, text='Permeability [mD]: '), 4, 0)
        self.layout().addWidget(self.permeability_text, 4, 1)
        self.layout().addWidget(QLabel(parent=self, text='krx: '), 5, 0)
        self.layout().addWidget(self.r_permeability_text, 5, 1)
        self.layout().addWidget(QLabel(parent=self, text='Resistance Factor Source: '), 6, 0)
        self.layout().addWidget(self.rf_combobox, 6, 1)
        self.offset_combobox.currentIndexChanged.connect(self.offset_source_select)
        self.permeability_combobox.currentIndexChanged.connect(self.permeability_source_select)
        self.rf_combobox.currentIndexChanged.connect(self.rf_source_select)

    def display_list_item_data(self):

        self.offset_combobox.setEnabled(True)
        self.permeability_combobox.setEnabled(True)
        self.rf_combobox.setEnabled(True)
        floods = self.parent().parent().parent().parent().experiment.binned_floods()
        index = self.flood_list.currentIndex()
        print(index)
        fi = floods[index.row()]

        for i, f in enumerate(floods):
            if fi.offset_source_flood == f:
                self.offset_combobox.blockSignals(True)
                self.offset_combobox.setCurrentIndex(i)
                self.offset_combobox.blockSignals(False)
            if fi.permeability_source_flood == f:
                self.permeability_combobox.blockSignals(True)
                self.permeability_combobox.setCurrentIndex(i)
                self.permeability_combobox.blockSignals(False)
            if fi.rf_source_flood == f:
                self.rf_combobox.blockSignals(True)
                self.rf_combobox.setCurrentIndex(i)
                self.rf_combobox.blockSignals(False)

        perm = fi.permeability
        ref_flood = self.parent().parent().parent().parent().experiment.petro_parameters['perm'][2]
        self.permeability_text.setText('{:.1f}, {:.1f}, {:.1f}, {:.1f}, {:.1f}'.format(*perm))
        if ref_flood is not None:
            ref_perm = ref_flood.permeability
            r_perm = divide(perm, ref_perm)
            self.r_permeability_text.setText('{:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.2f}'.format(*r_perm))

    def offset_source_select(self):

        floods = self.parent().parent().parent().parent().experiment.binned_floods()
        i = self.flood_list.currentIndex().row()
        j = self.offset_combobox.currentIndex()

        if floods[j].offset_source_flood == floods[i]:
            self.list_item_select()
            msg = 'Cannot have circular reference.'
            QMessageBox(parent=self, text=msg)
            return

        if floods[j].offset_source_flood != floods[j] and floods[j] != floods[i]:
            self.list_item_select()
            msg = 'Cannot set offset source to another flood with an offset source that is not itself.'
            QMessageBox(parent=self, text=msg)
            return

        floods[i].used_offsets = floods[j].used_offsets

        if floods[i].offset_source_flood is not floods[i]:
            floods[i].offset_source_flood.offset_referencer_floods.remove(floods[i])

        floods[i].offset_source_flood = floods[j]
        floods[j].offset_referencer_floods.append(floods[i])

        try:
            floods[i].reload_offsets_from_source()
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()

        cff = self.parent().parent().parent().parent().flood_forms[i]
        if cff is not None:
            print('Update core flood form offsets!')

    def permeability_source_select(self):

        floods = self.parent().parent().parent().parent().experiment.binned_floods()
        i = self.flood_list.currentIndex().row()
        j = self.permeability_combobox.currentIndex()

        if floods[j].permeability_source_flood == floods[i]:
            self.list_item_select()
            return

        if floods[j].permeability_source_flood != floods[j]:
            self.list_item_select()
            return

        floods[i].permeability = floods[j].permeability
        perm = floods[i].permeability
        ref_flood = self.parent().parent().parent().parent().experiment.petro_parameters['perm'][2]
        self.permeability_text.setText('{:.1f}, {:.1f}, {:.1f}, {:.1f}, {:.1f}'.format(*perm))
        if ref_flood is not None:
            ref_perm = ref_flood.permeability
            r_perm = divide(perm, ref_perm)
            self.r_permeability_text.setText('{:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.2f}'.format(*r_perm))

        if floods[i].permeability_source_flood is not floods[i]:
            floods[i].permeability_source_flood.permeability_referencer_floods.remove(floods[i])

        floods[i].permeability_source_flood = floods[j]
        floods[j].permeability_referencer_floods.append(floods[i])

    def rf_source_select(self):

        floods = self.parent().parent().parent().parent().experiment.binned_floods()
        i = self.flood_list.currentIndex().row()
        j = self.rf_combobox.currentIndex()

        try:
            if floods[j].rf_source_flood == floods[i]:
                self.list_item_select()
                return

            if floods[i].rf_source_flood is not floods[i]:
                floods[i].rf_source_flood.rf_referencer_floods.remove(floods[i])

            floods[i].rf_source_flood = floods[j]
            floods[j].rf_referencer_floods.append(floods[i])
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()


class EstimatesTab(QWidget):

    def __init__(self, parent):
        super(EstimatesTab, self).__init__(parent=parent)
        self.fe_view = parent.parent()

        lyt = QGridLayout()
        self.setLayout(lyt)
        font = self.font()
        font.setPointSize(11)
        self.setFont(font)
        txts = ['Init. State:', chr(954) + ' [md]:', chr(954) + 'rw0:', chr(954) + 'ro0:',
                'nw:', 'no:', 'Swr:', 'Sor:', 'SoI:', chr(966) + ':', 'C:']
        labels = []
        self.edits = []
        self.value_labels = []
        self.source_labels = []

        for i, txt in enumerate(txts):
            labels.append(QLabel(parent=self, text=txt))
            labels[-1].setAlignment(Qt.AlignRight)
            lyt.addWidget(labels[-1], i, 0, 1, 1)
            if i == 0:
                cb = QComboBox(parent=self)
                cb.addItems(['Brine', '2-phase', 'Oil'])
                cb.setCurrentIndex(self.fe_view.experiment.initial_state)
                cb.currentIndexChanged.connect(self.update_state)
                self.edits.append(cb)
            else:
                edt = QLineEdit(parent=self)
                edt.setFont(font)
                edt.setFixedWidth(100)
                self.value_labels.append(QLabel(parent=self, text='NaN'))
                self.value_labels[-1].setAlignment(Qt.AlignCenter)
                self.source_labels.append(QLabel(parent=self, text='None'))
                self.source_labels[-1].setAlignment(Qt.AlignCenter)
                edt.editingFinished.connect(self.update_values)
                self.edits.append(edt)

            lyt.addWidget(self.edits[-1], i, 1, 1, 1)

            if i > 0:
                lyt.addWidget(self.value_labels[-1], i, 2, 1, 1)
                lyt.addWidget(self.source_labels[-1], i, 3, 1, 1)
            else:
                v = QLabel(parent=self, text='Value:')
                v.setAlignment(Qt.AlignCenter)
                s = QLabel(parent=self, text='Source:')
                s.setAlignment(Qt.AlignCenter)
                lyt.addWidget(v, 0, 2, 1, 1)
                lyt.addWidget(s, 0, 3, 1, 1)

        self.update_state()
        self.set_values_and_sources()

    def set_values_and_sources(self):

        pps = self.fe_view.experiment.petro_parameters
        props = ['perm', 'krw0', 'kro0', 'nw', 'no', 'Swr', 'Sor', 'SoI', 'phi', 'C']

        for i, edit in enumerate(self.edits):

            edit.blockSignals(True)

            if i == 0:
                edit.setCurrentIndex(self.fe_view.experiment.initial_state)
            else:
                edit.setText(str(pps[props[i - 1]][0]))
                if props[i-1] == 'phi':
                    self.value_labels[i - 1].setText(str(round(pps[props[i - 1]][1], 4)))
                elif props[i-1] == 'C':
                    self.value_labels[i - 1].setText(str(round(pps[props[i - 1]][1], 1)))
                else:
                    self.value_labels[i - 1].setText(str(round(pps[props[i - 1]][1], 3)))
                source = pps[props[i - 1]][2]
                if source is None:
                    self.source_labels[i - 1].setText('Est.')
                else:
                    self.source_labels[i - 1].setText(source.name)

            edit.blockSignals(False)

    def update_values(self):

        ex = self.fe_view.experiment
        pps = ex.petro_parameters
        props = ['perm', 'krw0', 'kro0', 'nw', 'no', 'Swr', 'Sor', 'SoI', 'phi', 'C']

        for i, prop in enumerate(props):

            val0 = pps[prop][0]
            new_val = self.edits[i + 1].text()

            try:
                new_val = float(new_val)
            except ValueError as e:
                print(e)
                self.edits[i + 1].setText(str(val0))
                continue

            if pps[prop][0] == new_val:
                continue

            pps[prop][0] = new_val
            pps[prop][1] = new_val
            pps[prop][2] = None

        self.fe_view.experiment.petro_parameters = pps

    def update_state(self):

        self.fe_view.experiment.initial_state = self.edits[0].currentIndex()

        if self.fe_view.experiment.petro_parameters['SoI'][2] is None:
            if self.fe_view.experiment.initial_state == 0:
                self.fe_view.experiment.petro_parameters['SoI'][1] = 1.
            elif self.fe_view.experiment.initial_state == 2:
                self.fe_view.experiment.petro_parameters['SoI'][1] = 0.
            else:
                self.fe_view.experiment.petro_parameters['SoI'][1] = self.fe_view.experiment.petro_parameters['SoI'][0]

        # print(self.fe_view.experiment.petro_parameters['SoI'])
        self.fe_view.update_flood_phases()
        # print('setting values and sources.')
        self.set_values_and_sources()


class NotesTab(FloodListTab):

    def __init__(self, parent):
        super(NotesTab, self).__init__(parent=parent)

        self.notes_label = QLabel(parent=self, text='Planning Notes:')
        self.layout().addWidget(self.notes_label, 2, 0, 1, 2)
        self.notes_edit = QTextEdit(parent=self)
        self.layout().addWidget(self.notes_edit, 3, 0, 1, 2)
        self.observations_label = QLabel(parent=self, text='Observations:')
        self.layout().addWidget(self.observations_label, 4, 0, 1, 2)
        self.observations_edit = QTextEdit(parent=self)
        self.layout().addWidget(self.observations_edit, 5, 0, 1, 2)

        fl_list = parent.parent().parent().experiment.binned_floods()
        for fl in fl_list:
            self.flood_list.addItem(fl.name)

        self.list_item_selected.connect(self.display_list_item_data)
        self.notes_edit.textChanged.connect(lambda: self.set_list_item_data(0))
        self.observations_edit.textChanged.connect(lambda: self.set_list_item_data(1))

    def display_list_item_data(self):

        index = self.flood_list.currentIndex()
        floods = self.parent().parent().parent().parent().parent().experiment.binned_floods()
        fi = floods[index.row()]

        self.notes_edit.setText(fi.notes)
        self.observations_edit.setText(fi.observations)

    def set_list_item_data(self, notes_or_observations: int):

        index = self.flood_list.currentIndex()
        floods = self.parent().parent().parent().parent().parent().experiment.binned_floods()
        fi = floods[index.row()]

        if notes_or_observations == 0:
            fi.notes = self.notes_edit.toPlainText()
        else:
            fi.observations = self.observations_edit.toPlainText()


class ExperimentNotesTab(QWidget):

    def __init__(self, parent, sf: float):
        super(ExperimentNotesTab, self).__init__(parent=parent)

        self.objective_label = QLabel(parent=self, text='Objective Statement:')
        self.objective_edit = QTextEdit(parent=self)
        self.objective_edit.setFixedHeight(int(sf * 50))
        self.notes_label = QLabel('Planning Notes:')
        self.notes_edit = QTextEdit(parent=self)
        self.setLayout(QVBoxLayout())
        self.layout().addWidget(self.objective_label)
        self.layout().addWidget(self.objective_edit)
        self.layout().addWidget(self.notes_label)
        self.layout().addWidget(self.notes_edit)

        self.diplay_data()
        print('data displayed')
        self.objective_edit.textChanged.connect(lambda: self.set_data(0))
        self.notes_edit.textChanged.connect(lambda: self.set_data(1))
        print('experiment notes finished.')

    def diplay_data(self):

        expt = self.parent().parent().parent().parent().experiment
        self.notes_edit.setText(expt.notes)
        self.objective_edit.setText(expt.objective)

    def set_data(self, objective_or_notes: int):

        expt = self.parent().parent().parent().parent().parent().parent().experiment

        if objective_or_notes == 0:
            expt.objective = self.objective_edit.toPlainText()
        else:
            expt.notes = self.notes_edit.toPlainText()


class IFTTab(QWidget):

    def __init__(self, parent):
        super(IFTTab, self).__init__(parent=parent)

        lyt = QGridLayout()
        self.setLayout(lyt)

        self.oil_if_listbox = QListWidget(parent=self)
        self.brine_if_listbox = QListWidget(parent=self)
        ift_label = QLabel(parent=self, text='IFT [dynes/cm]:')
        self.ift_edit = QLineEdit(parent=self)
        self.ift_edit.setEnabled(False)
        self.ift_edit.editingFinished.connect(self.ift_edited)

        lyt.addWidget(self.oil_if_listbox, 0, 0, 1, 2)
        lyt.addWidget(self.brine_if_listbox, 1, 0, 1, 2)
        lyt.addWidget(ift_label, 2, 0, 1, 1)
        lyt.addWidget(self.ift_edit, 2, 1, 1, 1)

        self.oil_if_listbox.currentRowChanged.connect(self.check_enable_ift_edit)
        self.brine_if_listbox.currentRowChanged.connect(self.check_enable_ift_edit)

        self.oil_if_list = list()
        self.brine_if_list = list()
        self.experiment = self.parent().parent().experiment

        self.load_injection_fluids()

    def load_injection_fluids(self):

        self.oil_if_listbox.clear()
        self.brine_if_listbox.clear()

        expt = self.experiment

        oil_if_list = list()
        brine_if_list = list()

        self.oil_if_list = list()
        self.brine_if_list = list()

        for inj_fl in expt.injection_fluids_list.signal_lists[0].objects:
            if isinstance(inj_fl.specific_fluid, fluid.OilInjectionFluid):
                oil_if_list.append(inj_fl.name)
                self.oil_if_list.append(inj_fl)
                continue
            if isinstance(inj_fl.specific_fluid, fluid.BrineInjectionFluid):
                brine_if_list.append(inj_fl.name)
                self.brine_if_list.append(inj_fl)
                continue

        self.oil_if_listbox.addItems(oil_if_list)
        self.brine_if_listbox.addItems(brine_if_list)

    def check_enable_ift_edit(self, _):

        enable = False
        if self.oil_if_listbox.currentIndex().row() > -1 and self.brine_if_listbox.currentIndex().row() > -1:
            enable = True

        self.ift_edit.setEnabled(enable)

        if enable:
            self.load_ift()

    def load_ift(self):

        self.ift_edit.setText('')

        oi = self.oil_if_listbox.currentIndex().row()
        bi = self.brine_if_listbox.currentIndex().row()
        oif = self.oil_if_list[oi]
        bif = self.brine_if_list[bi]

        expt = self.experiment

        if (oif, bif) in expt.ift_dict.keys():
            self.ift_edit.setText('{:.3f}'.format(expt.ift_dict[(oif, bif)]))

    def set_ift(self, val: float):

        oi = self.oil_if_listbox.currentIndex().row()
        bi = self.brine_if_listbox.currentIndex().row()
        oif = self.oil_if_list[oi]
        bif = self.brine_if_list[bi]

        expt = self.experiment
        expt.ift_dict[(oif, bif)] = val

    def remove_ift(self):

        oi = self.oil_if_listbox.currentIndex().row()
        bi = self.brine_if_listbox.currentIndex().row()
        oif = self.oil_if_list[oi]
        bif = self.brine_if_list[bi]

        expt = self.experiment
        expt.ift_dict.pop((oif, bif))

    def ift_edited(self):

        txt = self.ift_edit.text()

        if not txt:
            self.remove_ift()
            return

        try:
            val = float(txt)
            if val < 0.:
                raise ValueError

            self.set_ift(val)

        except ValueError:

            self.load_ift()


class PetrophysicsView(QDialog):

    def __init__(self, experiment: FloodExperiment, parent: FloodExperimentView, sf: float=1.):
        super(PetrophysicsView, self).__init__(parent=parent)
        self.experiment = experiment
        self.setFixedSize(int(400 * sf), int(420 * sf))
        self.setWindowTitle('Petrophysics Tool')

        self.setLayout(QVBoxLayout())
        self.tab_widget = QTabWidget(parent=self)
        self.estimates_tab = EstimatesTab(self)
        self.permeability_tab = PermeabilityTab(self)
        self.ift_tab = IFTTab(self)
        # self.porosity_tab = PorosityTab(self)
        self.notes_tab_widget = QTabWidget(parent=self.tab_widget)
        self.experiment_notes_tab = ExperimentNotesTab(self.notes_tab_widget, sf)
        self.notes_tab = NotesTab(self.notes_tab_widget)
        self.tab_widget.addTab(self.estimates_tab, 'Est.')
        self.tab_widget.addTab(self.permeability_tab, chr(954))
        self.tab_widget.addTab(self.ift_tab, 'IFT')
        # self.tab_widget.addTab(self.porosity_tab, chr(966))
        # self.tab_widget.addTab(self.notes_tab, 'Notes')
        self.notes_tab_widget.addTab(self.experiment_notes_tab, 'Experiment')
        self.notes_tab_widget.addTab(self.notes_tab, 'Flood')
        self.tab_widget.addTab(self.notes_tab_widget, 'Notes')
        self.layout().addWidget(self.tab_widget)

        self.show()

    def closeEvent(self, a0: QCloseEvent):

        self.parent().petro_view = None


class SummaryViewer(QDialog):

    def __init__(self, parent: FloodExperimentView):
        super(SummaryViewer, self).__init__(parent=parent)
        self.sf = parent.sf
        self.experiment = parent.experiment

        self.resize(int(1100 * self.sf), int(600 * self.sf))
        # self.setFixedWidth(int(800 * self.sf))
        # self.setFixedHeight(int(400 * self.sf))
        self.setWindowTitle('Experiment Summary: {}'.format(self.experiment.name))

        main_lyt = QGridLayout()
        self.main_lyt = main_lyt
        self.setLayout(main_lyt)
        self.scroll_area = QScrollArea(parent=self)
        self.scroll_area.setWidgetResizable(True)
        self.scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.reload_button = QPushButton(parent=self, text='Reload')
        self.reload_button.clicked.connect(self.load_experiment)
        main_lyt.addWidget(self.reload_button)
        main_lyt.addWidget(self.scroll_area)

        lyt = QGridLayout()
        self.lyt = lyt
        self.scroll_widget = QWidget()
        self.scroll_area.setWidget(self.scroll_widget)
        self.scroll_widget.setLayout(lyt)

        self.core_properties_label = QLabel(parent=self.scroll_widget, text='Core Properties')
        self.core_name_label = QLabel(parent=self.scroll_widget)
        self.core_lithology_label = QLabel(parent=self.scroll_widget)
        self.core_diameter_label = QLabel(parent=self.scroll_widget)
        self.core_area_label = QLabel(parent=self.scroll_widget)
        self.core_length_label = QLabel(parent=self.scroll_widget)
        self.core_bulk_volume_label = QLabel(parent=self.scroll_widget)
        self.pore_volume_label = QLabel(parent=self.scroll_widget)
        self.porosity_label = QLabel(parent=self.scroll_widget)
        self.core_mass_label = QLabel(parent=self.scroll_widget)
        self.floods_label = QLabel(parent=self.scroll_widget, text='Flood Info')

        self.core_labels = [self.core_name_label, self.core_lithology_label, self.core_diameter_label,
                            self.core_area_label, self.core_length_label, self.core_bulk_volume_label,
                            self.pore_volume_label, self.porosity_label, self.core_mass_label]

        self.flood_labels = []
        self.plateau_labels = []

        font = self.core_name_label.font()
        font.setFamily('Courier')
        self.label_font = font
        font = self.core_name_label.font()
        font.setBold(True)
        font.setFamily('Courier')
        self.flood_label_font = font

        title_font = self.core_properties_label.font()
        title_font.setBold(True)
        title_font.setPointSize(12)
        self.title_font = title_font
        self.core_properties_label.setFont(title_font)
        self.core_properties_label.setAlignment(Qt.AlignHCenter)
        lyt.addWidget(self.core_properties_label)

        self.floods_label.setFont(title_font)
        self.floods_label.setAlignment(Qt.AlignHCenter)

        for i, label in enumerate(self.core_labels):
            label.setFont(font)
            lyt.addWidget(label, i + 1, 0, 1, 1)

        lyt.addWidget(self.floods_label)

        try:
            self.load_experiment()
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()

    def load_experiment(self):

        exp = self.experiment
        c = exp.core
        a = c.my_area()
        cl = c.my_length()
        cs = c.core_sections[0]
        d = sqrt((4. / pi) * a)
        porosity = exp.get_pv() / c.my_bulk_volume()

        self.core_name_label.setText('Core Name: {}'.format(c.name))
        self.core_lithology_label.setText('Lithology: {} ({} g/cc)'.format(cs.lithology, cs.densities[cs.lithology]))
        self.core_diameter_label.setText('D'.rjust(9) + ': {:.2f} cm ({:.3f} in)'.format(d, d / 2.54))
        self.core_area_label.setText('A'.rjust(9) + ': {:.2f} cm^2'.format(c.my_area()))
        self.core_length_label.setText('L'.rjust(9) + ': {:.2f} cm ({:.3f} in)'.format(cl, cl / 2.54))
        self.core_bulk_volume_label.setText('BV'.rjust(9) + ': {:.1f} cc'.format(c.my_bulk_volume()))
        self.pore_volume_label.setText('PV'.rjust(9) + ': {:.1f} mL'.format(exp.get_pv()))
        self.porosity_label.setText(chr(966).rjust(9) + ': {:.2f}%'.format(100. * porosity))
        self.core_mass_label.setText('M'.rjust(9) + ': {:.1f} g'.format(exp.get_core_mass(0)))

        for label in self.flood_labels:
            self.lyt.removeWidget(label)
            label.close()

        for label in self.plateau_labels:
            self.lyt.removeWidget(label)
            label.close()

        perm_ref_flood = self.experiment.petro_parameters['perm'][2]
        temp_perm = self.experiment.petro_parameters['perm'][0]
        ref_perm = array([temp_perm for ii in range(5)])
        if perm_ref_flood is not None:
            ref_perm = perm_ref_flood.permeability

        sw0 = exp.get_initial_saturation()
        rho_o = np.nan
        rho_w = np.nan
        sw_final = 1.
        row = 11

        for i, fl in enumerate(exp.binned_floods()):
            label = QLabel(parent=self.scroll_widget)
            if i == 0:
                sw_init = sw0
                sw_final = exp.get_flood_saturation(fl)
            else:
                sw_init = sw_final
                sw_final = exp.get_flood_saturation(fl)

            flood_text = '{}, T = {:.1f}' + chr(176) + 'C, Sw {:.3f} -> {:.3f}'
            label.setText(flood_text.format(fl.name, fl.temperature, sw_init, sw_final))
            if isinstance(fl.fluid.specific_fluid, fluid.OilInjectionFluid) or fl.get_background_color() == 'black':
                label.setStyleSheet('QLabel{color: white; background-color: ' + fl.get_background_color() + '};')
            else:
                label.setStyleSheet('QLabel{background-color: ' + fl.get_background_color() + '};')
            label.setFont(self.flood_label_font)
            self.lyt.addWidget(label, row, 0, 1, 1)
            self.flood_labels.append(label)

            row += 1

            p_flow_rates, _ = fl.get_rates_and_plateaus(0)
            p_fluids = fl.get_plateau_fluids([])

            if fl.plateau_whole:

                use_rf = False
                if len(fl.plateau_whole_rf) == len(fl.plateau_whole):
                    use_rf = True

                j = 0
                for plateau_w, plateau_1, plateau_2, plateau_3, plateau_4 in zip(fl.plateau_whole, fl.plateau_sec1,
                                                                                 fl.plateau_sec2, fl.plateau_sec3,
                                                                                 fl.plateau_sec4):
                    p_label = QLabel(parent=self.scroll_widget)
                    p_label.setFont(self.label_font)
                    q = p_flow_rates[j]
                    mu = fl.get_fluid_viscosity(j)
                    p_text = '    ' + \
                             'W: ' + '{:.1f}'.format(plateau_w).rjust(5) + ' psi' \
                             '  1: ' + '{:.1f}'.format(plateau_1).rjust(5) + ' psi' \
                             '  2: ' + '{:.1f}'.format(plateau_2).rjust(5) + ' psi' \
                             '  3: ' + '{:.1f}'.format(plateau_3).rjust(5) + ' psi' \
                             '  4: ' + '{:.1f}'.format(plateau_4).rjust(5) + ' psi'

                    p_label.setText(p_text)
                    p_label.setStyleSheet('QLabel{background-color: yellow};')
                    self.lyt.addWidget(p_label, row, 0, 1, 1)
                    self.plateau_labels.append(p_label)

                    p_label_right = QLabel(parent=self.scroll_widget)
                    p_label_right.setStyleSheet('QLabel{background-color: yellow};')
                    p_label_right.setFont(self.flood_label_font)
                    p_text_right = '  Q:  {:.3f} mL/min  |  '.format(q) + chr(956) + ':' + \
                                   '{:.2f} cP'.format(mu).rjust(11)
                    p_label_right.setText(p_text_right)
                    self.lyt.addWidget(p_label_right, row, 1, 1, 1)
                    self.plateau_labels.append(p_label_right)

                    row += 1

                    dpis = [plateau_w, plateau_1, plateau_2, plateau_3, plateau_4]
                    clis = [cl, 0, 0, 0, 0]
                    for k in range(4):
                        if cl >= 7.62 * k:
                            clis[k + 1] = 7.62
                            continue
                        clis[k + 1] = cl - sum(clis[1:k+1])
                        break

                    p_label = QLabel(parent=self.scroll_widget)
                    p_label.setFont(self.flood_label_font)
                    perm = [(245. * q * mu * cli) / (a * dpi) for cli, dpi in zip(clis, dpis)]
                    ftperd = q * 60. * 24. * 12. / exp.get_pv() / (cl / 2.54)
                    p_text = '    ' + \
                             ' ' + '{:.1f}'.format(perm[0]).rjust(7) + ' mD ' \
                             '   ' + '{:.1f}'.format(perm[1]).rjust(7) + ' mD ' \
                             '   ' + '{:.1f}'.format(perm[2]).rjust(7) + ' mD ' \
                             '   ' + '{:.1f}'.format(perm[3]).rjust(7) + ' mD ' \
                             '   ' + '{:.1f}'.format(perm[4]).rjust(7) + ' mD '

                    p_label.setText(p_text)
                    p_label.setStyleSheet('QLabel{background-color: yellow; color: red};')
                    self.lyt.addWidget(p_label, row, 0, 1, 1)
                    self.plateau_labels.append(p_label)

                    p_label_right = QLabel(parent=self.scroll_widget)
                    p_label_right.setStyleSheet('QLabel{background-color: yellow};')
                    p_label_right.setFont(self.flood_label_font)
                    p_text_right = '      {:.2f} ft/d  '.format(ftperd)
                    p_label_right.setText(p_text_right)
                    self.lyt.addWidget(p_label_right, row, 1, 1, 1)
                    self.plateau_labels.append(p_label_right)

                    row += 1

                    p_label = QLabel(parent=self.scroll_widget)
                    p_label.setFont(self.flood_label_font)
                    if isinstance(fl.fluid.specific_fluid, fluid.OilInjectionFluid):
                        p_text = '  kro:  '
                    elif fl.fluid.get_fluid_type() == fluid.FluidType.POLYMER_SOLUTION:
                        p_text = '  krp:  '
                    else:
                        p_text = '  krw:  '

                    for permi, ref_permi in zip(perm, ref_perm):
                        p_text += '{:.2f}'.format(permi / ref_permi).ljust(14)

                    p_label.setText(p_text)
                    p_label.setStyleSheet('QLabel{background-color: yellow; color: red};')
                    self.lyt.addWidget(p_label, row, 0, 1, 1)
                    self.plateau_labels.append(p_label)

                    if isinstance(fl.fluid.specific_fluid, fluid.OilInjectionFluid):
                        rho_o = fl.fluid.specific_fluid.get_density(fl.temperature)

                    if isinstance(fl.fluid.specific_fluid, fluid.BrineInjectionFluid):
                        rho_w = fl.fluid.specific_fluid.get_density(fl.temperature)

                    if fl.n_phases < 2 or isinstance(fl.fluid.specific_fluid, fluid.OilInjectionFluid):
                        row += 1
                        j += 1
                        continue

                    ift = exp.get_ift(fl, p_fluids[j])

                    nc = mu * (ftperd * porosity * 12. * 2.54 / 24. / 60. / 60.) / ift

                    nb = -980 * (rho_o - rho_w) * ref_perm[0] * 9.869e-12 / ift

                    trapping_n = nc
                    if not np.isnan(nb):
                        if not fl.reverse_flow:
                            trapping_n += nb
                        else:
                            trapping_n -= nb

                    p_label_right = QLabel(parent=self.scroll_widget)
                    p_label_right.setStyleSheet('QLabel{background-color: yellow};')
                    p_label_right.setFont(self.flood_label_font)
                    p_text_right = '  NT: {:.2E} (IFT = {:.1E} dynes/cm)'.format(trapping_n, ift)
                    p_label_right.setText(p_text_right)
                    self.lyt.addWidget(p_label_right, row, 1, 1, 1)
                    self.plateau_labels.append(p_label_right)

                    row += 1

                    if not use_rf:
                        j += 1
                        continue

                    p_label = QLabel(parent=self.scroll_widget)
                    p_label.setFont(self.flood_label_font)
                    p_text = '   RF:  '
                    rfs = [fl.plateau_whole_rf[j], fl.plateau_sec1_rf[j], fl.plateau_sec2_rf[j],
                           fl.plateau_sec3_rf[j], fl.plateau_sec4_rf[j]]

                    for rf in rfs:
                        p_text += '{:.2f}'.format(rf).ljust(14)

                    p_label.setText(p_text)
                    p_label.setStyleSheet('QLabel{background-color: yellow; color: red};')
                    self.lyt.addWidget(p_label, row, 0, 1, 1)
                    self.plateau_labels.append(p_label)

                    row += 1

                    j += 1

            label = QLabel(parent=self.scroll_widget)
            label.setFont(self.flood_label_font)
            self.lyt.addWidget(label, row, 0, 1, 1)
            self.flood_labels.append(label)

            row += 1

    def closeEvent(self, a0: QCloseEvent):

        self.parent().summary_viewer = None


if __name__ == "__main__":
    import sys
    from PyQt5.Qt import QApplication

    app = QApplication(sys.argv)
    core = core.Core('', [3.76], [5.88], [137.99], 'sandstone')
    experiment = FloodExperiment('', core)
    view = FloodExperimentView(experiment)

    sys.exit(app.exec_())
