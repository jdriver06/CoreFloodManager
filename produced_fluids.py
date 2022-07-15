
from PyQt5.QtWidgets import QWidget, QTableWidget, QTableWidgetItem, QHBoxLayout, QVBoxLayout, QLabel, \
    QLineEdit, QDialog, QPushButton, QGridLayout, QComboBox, QCheckBox, QMessageBox, QRadioButton, QButtonGroup, \
    QFrame, QGraphicsEllipseItem
from PyQt5.Qt import QIcon, pyqtSignal, Qt, QFont, QApplication, QPoint, QPointF
from PyQt5.QtGui import QCloseEvent, QKeyEvent, QMouseEvent
from PyQt5.QtChart import QChart, QChartView, QValueAxis, QLineSeries, QScatterSeries
from PyQt5.QtCore import QItemSelection
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import j_utils as utils
from copy import copy
from core import CoreSection
from fluid import FluidType, Brine
import flood


__version__ = '0.1.0'


class ProducedFluids:

    def __init__(self, fld, aq_phase: bool, oil_phase: bool):

        if not aq_phase and not oil_phase:
            raise Exception

        self.flood = fld

        if aq_phase and oil_phase:
            self.volumes = np.array([[], [], [], []])
        else:
            self.volumes = np.array([[]])

        self.revisions = np.array([[], [], [], []])
        self.micro_emulsions = np.array([[], []])
        self.measurements = np.array([[], [], [], [], []])
        self.ions = np.array([[], [], [], [], [], [], [], [], [], []])
        self.ppm_to_cp = np.array([[], []])
        self.ppm_to_abs = np.array([[], []])
        self.ppm_to_cp_fit_params = np.array([])
        self.ppm_to_abs_fit_params = np.array([])
        self.cp_tracer_object = None
        self.abs_tracer_object = None
        self.cp_normalized = True
        self.abs_normalized = True
        self.cp_or_abs_used = 0
        self.ppm = np.array([])
        self.retention = np.nan
        self.retention_ppm_min = 0.

        self.aq_phase = aq_phase
        self.oil_phase = oil_phase
        self.accumulation_time = 60.
        self.volume_factor = 1.
        self.swelling_factor = 1.
        self.gamma = 7.3
        self.temperature = 23.
        self.read_temperature = 23.
        self.dead_volume = 0.

    def micro_emulsion_oil_or_water(self, is_oil: bool) -> np.array:

        me = self.micro_emulsions
        same = me.shape[1] == self.volumes.shape[1]
        if same:
            me_type = me[0, :]
            me_frac = me[1, :]
        else:
            me_type = me[0, :-1]
            me_frac = me[1, :-1]
        if not is_oil:
            me_frac = 1. - me_frac

        same = self.volumes.shape[1] == self.revisions.shape[1]
        if same:
            v = self.volumes + self.revisions
        else:
            v = self.volumes + self.revisions[:, :-1]
        water = v[0, :]
        oil = v[3, :]
        t3me = v[1, :] - v[0, :]

        me_oil_or_water = np.zeros((v.shape[1], ))

        type_i = np.where(me_type == 1.)
        type_ii = np.where(me_type == 2.)
        type_iii = np.where(me_type == 3.)

        if type_i[0].size:
            # print('There are Type I micro-emulsions.')
            me_oil_or_water[type_i] = np.multiply(me_frac[type_i], water[type_i])
        if type_ii[0].size:
            # print('There are Type II micro-emulsions.')
            me_oil_or_water[type_ii] = np.multiply(me_frac[type_ii], oil[type_ii])
        if type_iii[0].size:
            # print('There are Type III micro-emulsions.')
            me_oil_or_water[type_iii] = np.multiply(me_frac[type_iii], t3me[type_iii])

        return me_oil_or_water

    def micro_emulsion_oil(self) -> np.array:

        return self.micro_emulsion_oil_or_water(True)

    def micro_emulsion_water(self) -> np.array:

        return self.micro_emulsion_oil_or_water(False)

    def thermal_correction_water(self) -> float:

        di = Brine('DI')
        rho_read = di.get_density(self.read_temperature)
        rho_flood = di.get_density(self.flood.temperature)

        print(rho_read / rho_flood)

        return rho_read / rho_flood

    def thermal_correction_oil(self):

        oil_if = self.flood.experiment.get_most_recent_oil_if(self.flood)
        if oil_if is None:
            return 1.

        oil = oil_if.specific_fluid.oil_sample.ref_objects[0]
        ft = oil_if.get_fluid_type()

        if ft == FluidType.OIL:
            rho_read = oil.get_density(self.read_temperature)
            rho_flood = oil.get_density(self.flood.temperature)
            return rho_read / rho_flood

        elif ft == FluidType.LIVE_OIL:
            rho_read = oil.get_density(self.read_temperature)
            rho_flood = oil.get_live_oil_density(self.flood.back_pressure)
            return rho_read / rho_flood

        return 1.

    def final_oil_and_water(self) -> np.array:

        if self.volumes.shape[1] == 0:
            return np.array([[], [], []])

        final = np.zeros((3, self.volumes.shape[1]))

        vf = self.volume_factor
        sf = self.swelling_factor

        same = self.revisions.shape[1] == self.volumes.shape[1]
        if same:
            regular_oil = self.volumes[3, :].copy() + self.revisions[3, :].copy()
            type_ii = np.where(self.micro_emulsions[0, :] == 2.)
        else:
            regular_oil = self.volumes[3, :].copy() + self.revisions[3, :-1].copy()
            type_ii = np.where(self.micro_emulsions[0, :-1] == 2.)
        if type_ii[0].size != 0:
            regular_oil[type_ii] = 0.

        oil = regular_oil + self.micro_emulsion_oil()
        water_in_oil = oil * (1. - 1. / sf)
        oil -= water_in_oil
        oil *= vf

        print('oil_if thermal correction = {:.3f}'.format(self.thermal_correction_oil()))

        if same:
            regular_water = self.volumes[0, :].copy() + self.revisions[0, :].copy()
            type_i = np.where(self.micro_emulsions[0, :] == 1.)
        else:
            regular_water = self.volumes[0, :].copy() + self.revisions[0, :-1].copy()
            type_i = np.where(self.micro_emulsions[0, :-1] == 1.)
        if type_i[0].size != 0:
            regular_water[type_i] = 0.

        water = regular_water + water_in_oil + self.micro_emulsion_water()
        water *= self.thermal_correction_water()

        final[0, :] = water
        final[2, :] = oil
        final[1, :] = water + oil

        # print('final:', final)

        return final

    def oil_cut(self):

        final_ow = self.final_oil_and_water()
        return np.divide(final_ow[2, :], final_ow[1, :])

    def cumulative_volume(self):

        final_ow = self.final_oil_and_water()
        return np.cumsum(final_ow[1, :])

    def cumulative_volume_centered(self):

        final_ow = self.final_oil_and_water()
        return np.cumsum(final_ow[1, :]) - 0.5 * final_ow[1, :]

    def cumulative_oil_volume(self):

        final_ow = self.final_oil_and_water()
        return np.cumsum(final_ow[2, :])

    def cumulative_aq_volume(self):

        final_ow = self.final_oil_and_water()
        return np.cumsum(final_ow[0, :])

    def oil_produced(self):

        final_ow = self.final_oil_and_water()
        if final_ow.shape[1] == 0:
            return 0.
        return np.sum(final_ow[2, :])

    def brine_produced(self):
        final_ow = self.final_oil_and_water()
        if final_ow.shape[1] == 0:
            return 0.
        return np.sum(final_ow[0, :])

    def plot_oil_cut(self):

        oc = self.oil_cut()
        cv = self.cumulative_volume_centered()

        plt.plot(cv, oc)
        plt.title(str(np.round(self.oil_produced(), 3)) + ' mL oil produced.')
        plt.xlabel('Injected Volume [mL]')
        plt.show()


class ProducedFluidsView(QDialog):

    def __init__(self, pf: ProducedFluids, parent=None):

        if not pf.aq_phase and not pf.oil_phase:
            return

        super(ProducedFluidsView, self).__init__(parent=parent)

        if parent is None:
            self.setWindowTitle('Produced Fluids Tool')
        else:
            self.setWindowTitle('Produced Fluids Tool: ' + parent.flood.name)

        self.setWindowIcon(QIcon('Coreholder Cropped.jpg'))
        if parent is not None:
            sf = parent.width() / 800.
        else:
            sf = 1.
        self.resize(int(sf * 1080), int(sf * 500))
        self.setFixedWidth(int(sf * 1080))
        layout = QVBoxLayout()
        self.setLayout(layout)

        font = self.font()
        font.setPointSize(10)

        widget = QWidget(parent=self)
        # layout.addWidget(widget)
        w_layout = QGridLayout()
        widget.setLayout(w_layout)
        # text = QLabel(parent=self, text='Accumulation Time [min]:')
        # text.setFont(font)
        # w_layout.addWidget(text)
        # self.accumulation_time_edit = QLineEdit(parent=self)
        # self.accumulation_time_edit.setText(str(pf.accumulation_time))
        # self.accumulation_time_edit.setFont(font)
        # self.accumulation_time_edit.editingFinished.connect(self.time_update)
        # w_layout.addWidget(self.accumulation_time_edit)

        text = QLabel(parent=self, text='Read Temperature:')
        text.setFont(font)
        w_layout.addWidget(text, 0, 0, 1, 1)
        self.read_temperature_edit = QLineEdit(parent=self)
        self.read_temperature_edit.setText(str(pf.read_temperature))
        self.read_temperature_edit.setFont(font)
        self.read_temperature_edit.editingFinished.connect(self.rt_update)
        w_layout.addWidget(self.read_temperature_edit, 0, 1, 1, 1)
        text = QLabel(parent=self, text='Volume Factor:')
        text.setFont(font)
        w_layout.addWidget(text, 1, 0, 1, 1)
        self.volume_factor_edit = QLineEdit(parent=self)
        self.volume_factor_edit.setText(str(pf.volume_factor))
        self.volume_factor_edit.setFont(font)
        self.volume_factor_edit.editingFinished.connect(self.vf_update)
        w_layout.addWidget(self.volume_factor_edit, 1, 1, 1, 1)
        text = QLabel(parent=self, text='Swelling Factor:')
        text.setFont(font)
        w_layout.addWidget(text, 2, 0, 1, 1)
        self.swelling_factor_edit = QLineEdit(parent=self)
        self.swelling_factor_edit.setText(str(pf.swelling_factor))
        self.swelling_factor_edit.setFont(font)
        self.swelling_factor_edit.editingFinished.connect(self.sf_update)
        w_layout.addWidget(self.swelling_factor_edit, 2, 1, 1, 1)
        self.suggest_vf_button = QPushButton(parent=self, text='Suggest VF')
        self.suggest_vf_button.setFont(font)
        self.suggest_vf_button.clicked.connect(self.suggest_vf)
        w_layout.addWidget(self.suggest_vf_button, 1, 2, 1, 1)
        self.oil_cut_button = QPushButton(parent=self, text='Oil Cut')
        self.oil_cut_button.setFont(font)
        self.oil_cut_button.clicked.connect(self.show_oil_cut)
        w_layout.addWidget(self.oil_cut_button, 2, 2, 1, 1)

        m_widget = QWidget(parent=self)
        m_lyt = QGridLayout()
        m_widget.setLayout(m_lyt)
        text = QLabel(parent=self, text=u' Shear rate for ' + chr(951) + u' [1/s]:')
        text.setFont(font)
        m_lyt.addWidget(text, 0, 0, 1, 1)
        self.shear_rate_edit = QLineEdit(parent=self)
        self.shear_rate_edit.setText(str(pf.gamma))
        self.shear_rate_edit.setFont(font)
        self.shear_rate_edit.editingFinished.connect(self.gamma_update)
        m_lyt.addWidget(self.shear_rate_edit, 0, 1, 1, 1)
        text = QLabel(parent=self, text=u' Temp. for ' + chr(951) + u' [' + chr(176) + 'C]:')
        text.setFont(font)
        m_lyt.addWidget(text, 0, 2, 1, 1)
        self.temp_edit = QLineEdit(parent=self)
        self.temp_edit.setText(str(pf.temperature))
        self.temp_edit.setFont(font)
        self.temp_edit.editingFinished.connect(self.temp_update)
        m_lyt.addWidget(self.temp_edit, 0, 3, 1, 1)
        self.porosity_button = QPushButton(parent=self, text='Porosity')
        self.porosity_button.setFont(font)
        self.porosity_button.clicked.connect(self.porosity_tool)
        m_lyt.addWidget(self.porosity_button, 0, 4, 1, 1)
        self.retention_button = QPushButton(parent=self, text='Retention')
        self.retention_button.setFont(font)
        self.retention_button.clicked.connect(self.retention_tool)
        m_lyt.addWidget(self.retention_button, 0, 5, 1, 1)
        self.alkali_consumption_button = QPushButton(parent=self, text='Alkali Consumption')
        self.alkali_consumption_button.setFont(font)
        self.alkali_consumption_button.clicked.connect(self.alkali_consumption_tool)
        m_lyt.addWidget(self.alkali_consumption_button, 1, 4, 1, 2)

        g_layout = QGridLayout()
        layout.addLayout(g_layout)

        self.produced_fluids = pf

        self.volumes_radio_group = QButtonGroup(parent=self)
        self.volumes_frame = QFrame(parent=self)
        self.volumes_frame.setLayout(QHBoxLayout())
        v_names = ['Original', 'Revisions', 'Net', 'Micro-emulsions', 'Final']
        self.volumes_radios = []
        for name in v_names:
            radio = QRadioButton(parent=self, text=name)
            self.volumes_radios.append(radio)
            self.volumes_radio_group.addButton(radio)
            self.volumes_frame.layout().addWidget(radio)

        self.measurements_radio_group = QButtonGroup(parent=self)
        self.measurements_frame = QFrame(parent=self)
        self.measurements_frame.setLayout(QHBoxLayout())
        m_names = ['Standard', 'Ions', 'Surfactants']
        self.measurements_radios = []
        for name in m_names:
            radio = QRadioButton(parent=self, text=name)
            self.measurements_radios.append(radio)
            self.measurements_radio_group.addButton(radio)
            self.measurements_frame.layout().addWidget(radio)

        self.table = ProducedFluidsTable(parent=self)
        self.revisions_table = RevisionsTable(parent=self)
        self.revisions_table.setVisible(False)
        self.net_table = NetTable(parent=self)
        self.net_table.setVisible(False)
        self.micro_emulsions_table = METable(parent=self)
        self.micro_emulsions_table.setVisible(False)
        self.final_values_table = FinalValuesTable(parent=self)
        self.final_values_table.setVisible(False)
        self.measurements_table = EffluentMeasurementsTable(parent=self)
        self.ions_table = IonsTable(parent=self)
        self.ions_table.setVisible(False)
        self.table.rc_changed.connect(self.revisions_table.sync_rows)
        self.table.data_changed.connect(self.revisions_table.reset_data_at_row)
        self.table.rc_changed.connect(self.net_table.sync_rows)
        self.table.rc_changed.connect(self.micro_emulsions_table.sync_rows)
        self.table.rc_changed.connect(self.final_values_table.sync_rows)
        self.revisions_table.data_changed.connect(self.micro_emulsions_table.update_data_at_row)
        self.table.rc_changed.connect(self.measurements_table.sync_rows)
        self.table.rc_changed.emit(self.table.rowCount())

        g_layout.addWidget(self.volumes_frame, 0, 0)
        g_layout.addWidget(self.measurements_frame, 0, 1)
        g_layout.addWidget(self.table, 1, 0)
        g_layout.addWidget(self.measurements_table, 1, 1)
        # widget.setFixedWidth(self.table.width())
        g_layout.addWidget(widget, 2, 0)
        m_widget.setFixedWidth(self.measurements_table.width())
        g_layout.addWidget(m_widget, 2, 1)
        print(widget.width(), self.table.width(), self.volumes_frame.width(), m_widget.width())

        self.g_layout = g_layout

        self.volumes_radios[0].setChecked(True)
        self.measurements_radios[0].setChecked(True)

        self.volumes_radio_group.buttonToggled.connect(self.volumes_radio_toggled)
        self.measurements_radio_group.buttonToggled.connect(self.measurements_radio_toggled)
        self.measurements_radios[-1].setEnabled(False)

        self.r_tool = None
        self.p_tool = None
        self.c_tool = None
        self.tracer_table = self.measurements_table
        self.tracer_target = np.nan

        if isinstance(parent.flood, flood.MultiCoreFlood):
            self.porosity_button.setEnabled(False)

        self.show()

    def radio_toggled(self, tables: list, ind: int, col: int):

        visibles = []
        for table in tables:
            visibles.append(table.isVisible())
        v_ind = visibles.index(True)

        if ind == v_ind:
            return

        i = -1
        for table, visible in zip(tables, visibles):
            i += 1
            if visible or i != ind:
                table.setVisible(False)
                self.g_layout.removeWidget(table)
            else:
                table.setVisible(True)
                self.g_layout.addWidget(table, 1, col)

    @staticmethod
    def find_checked(radios: list) -> QRadioButton:

        checked = None
        for radio in radios:
            if radio.isChecked():
                checked = radio
                break

        return checked

    def volumes_radio_toggled(self):

        checked = self.find_checked(self.volumes_radios)
        if checked is None:
            return

        names = ['Original', 'Revisions', 'Net', 'Micro-emulsions', 'Final']
        ind = names.index(checked.text())
        tables = [self.table, self.revisions_table, self.net_table, self.micro_emulsions_table, self.final_values_table]

        self.radio_toggled(tables, ind, 0)

    def measurements_radio_toggled(self):

        checked = self.find_checked(self.measurements_radios)
        if checked is None:
            return

        names = ['Standard', 'Ions']
        ind = names.index(checked.text())
        tables = [self.measurements_table, self.ions_table]

        self.radio_toggled(tables, ind, 1)
        self.tracer_table = tables[ind]
        if ind == 1:
            self.tracer_target = 1
        else:
            self.tracer_target = np.nan

    def suggest_vf(self):

        vf = self.produced_fluids.thermal_correction_oil()
        if np.isnan(vf):
            vf = 1.

        self.volume_factor_edit.setText(str(np.round(vf, 3)))
        self.vf_update()

    def show_oil_cut(self):

        s = self.produced_fluids.volumes.shape
        if not s[1]:
            return

        cvc = self.produced_fluids.cumulative_volume_centered()
        oc = self.produced_fluids.oil_cut()
        fld = self.parent().flood
        expt = fld.experiment
        flds_list = expt.binned_floods(fld)
        ind = flds_list.index(fld)

        if not ind:
            return

        sw = expt.get_flood_saturation(flds_list[ind - 1])
        soi = 1. - sw
        pv = self.parent().flood.get_experiment_total_pore_volume()
        so = soi * np.ones(cvc.shape)
        cov = self.produced_fluids.cumulative_oil_volume()
        cor = cov / (pv * soi)
        so = so - cov / pv
        cvc2 = np.zeros((cvc.shape[0] + 1, ))
        cvc2[1:] = cvc[:]
        so2 = np.zeros(cvc2.shape)
        so2[0] = soi
        so2[1:] = so[:]
        cor2 = np.zeros(cvc2.shape)
        cor2[1:] = cor[:]

        title_font = {'fontname': 'Arial'}
        font = font_manager.FontProperties(family='Arial')
        axis_font = {'font.family': 'monospace', 'font.monospace': 'Courier New'}

        plt.rcParams.update(axis_font)

        plt.plot(cvc / pv, 100. * oc, 'r')
        plt.plot(cvc2 / pv, 100. * so2, 'g')
        plt.plot(cvc2 / pv, 100. * cor2, 'b')

        title_str = '{:.3f} mL oil produced, Sor = {:.1f}%, COR = {:.1f}%'.format(cov[-1],
                                                                                  100. * so[-1],
                                                                                  100. * cor[-1])
        plt.title(title_str, **title_font)
        plt.xlabel('Pore Volumes Injected', **title_font)
        plt.ylabel('[%]', **title_font)
        plt.legend(['Oil Cut', 'So', 'Cum. Oil Rec.'], prop=font)
        plt.ylim(0., 100.)
        _, x_max = plt.gca().get_xlim()
        plt.xlim(0., x_max)
        plt.show()

    def porosity_tool(self):
        print(self.p_tool)
        if self.p_tool is not None:
            return
        try:
            self.p_tool = PorosityTool(self)
            if not self.p_tool.isVisible():
                self.p_tool.close()
        except Exception as e:
            print(e)

    def retention_tool(self):
        try:
            if self.r_tool is None:
                self.r_tool = RetentionTool(self)
        except Exception as e:
            print(e)
            # msg = QMessageBox(parent=self, text=e)
            # msg.exec_()

    def alkali_consumption_tool(self):
        try:
            if self.c_tool is None:
                self.c_tool = AlkaliConsumptionTool(self)

        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()

    def time_update(self):

        time = self.value_update(self.accumulation_time_edit)

        if time == -1:
            self.accumulation_time_edit.setText(str(self.produced_fluids.accumulation_time))
        else:
            self.produced_fluids.accumulation_time = time

        # self.show_values()

    def rt_update(self):

        rt = self.value_update(self.read_temperature_edit)

        if rt == -1:
            self.read_temperature_edit.setText(str(self.produced_fluids.read_temperature))
        else:
            self.produced_fluids.read_temperature = rt

    def vf_update(self):

        vf = self.value_update(self.volume_factor_edit)

        if vf == -1:
            self.volume_factor_edit.setText(str(self.produced_fluids.volume_factor))
        else:
            self.produced_fluids.volume_factor = vf

        # self.show_values()

    def sf_update(self):

        sf = self.value_update(self.swelling_factor_edit)

        if sf == -1:
            self.swelling_factor_edit.setText(str(self.produced_fluids.swelling_factor))
        else:
            self.produced_fluids.swelling_factor = sf

    def gamma_update(self):

        gamma = self.value_update(self.shear_rate_edit)

        if gamma == -1:
            self.shear_rate_edit.setText(str(self.produced_fluids.gamma))
        else:
            self.produced_fluids.gamma = gamma

        # self.show_values()

    def temp_update(self):

        temp = self.value_update(self.temp_edit)

        if temp == -1:
            self.temp_edit.setText(str(self.produced_fluids.temperature))
        else:
            self.produced_fluids.temperature = temp

        # self.show_values()

    @staticmethod
    def value_update(source_edit: QLineEdit):

        value = source_edit.text()

        try:
            value = float(value)
            if value <= 0:
                raise ValueError

        except ValueError:
            return -1

        return value

    def show_values(self):

        pf = self.produced_fluids
        print('values: ', pf.accumulation_time, pf.volume_factor, pf.gamma, pf.temperature)

    def closeEvent(self, a0: QCloseEvent):
        if self.r_tool is not None:
            self.r_tool.close()
        if self.p_tool is not None:
            self.p_tool.close()
        if self.parent() is not None:
            self.parent().effluent_view = None
            print(self.produced_fluids.measurements[1, :])
            self.parent().flood.effluent = self.produced_fluids
            self.parent().linked_widget.update_flood_forms_saturation_text()


class ProducedFluidsTable(QTableWidget):

    rc_changed = pyqtSignal(int)
    data_changed = pyqtSignal(int)

    def __init__(self, parent: ProducedFluidsView):
        super(ProducedFluidsTable, self).__init__(parent=parent)

        self.control_pressed = False

        font = self.font()
        font.setFamily('Courier')
        font.setPointSize(12)
        self.setFont(font)

        self.setColumnCount(4)
        self.setHorizontalHeaderLabels(['Aq [mL]', 'M.I. [mL]', 'Tot [mL]', 'Oil [mL]'])

        self.sync_to_data_vols()

        self.itemChanged.connect(self.check_item)

    def sync_to_data_vols(self):

        self.sync_rows()
        self.sync_to_data(self.parent().produced_fluids.volumes)

    def sync_rows(self):

        rc_old = self.rowCount()
        volumes = self.parent().produced_fluids.volumes
        rc_new = volumes.shape[1] + 1

        if rc_new == rc_old:
            return

        self.setRowCount(rc_new)

        if rc_new < rc_old:
            for col in range(self.columnCount()):
                if self.item(rc_new - 1, col) is not None:
                    item = self.takeItem(rc_new - 1, col)
                    del item

        self.rc_changed.emit(rc_new)

    def sync_to_data(self, data: np.array):

        for r in range(data.shape[1]):
            for c in range(data.shape[0]):

                datum = data[c, r]
                txt = ''
                if not np.isnan(datum):
                    txt = str(np.round(datum, 3))

                item = self.item(r, c)
                self.blockSignals(True)
                if item is None:
                    item = QTableWidgetItem('')
                    self.setItem(r, c, item)
                item.setText(txt)
                self.blockSignals(False)

    def sync_oil(self, r: int):

        volumes = self.parent().produced_fluids.volumes
        if np.isnan(volumes[1, r]):
            volumes[3, r] = np.round(volumes[2, r] - volumes[0, r], 3)
        else:
            volumes[3, r] = np.round(volumes[2, r] - volumes[1, r], 3)

    def sync_total(self, r: int):

        volumes = self.parent().produced_fluids.volumes
        if np.isnan(volumes[1, r]):
            volumes[2, r] = np.round(volumes[0, r] + volumes[3, r], 3)
        else:
            volumes[2, r] = np.round(volumes[1, r] + volumes[3, r], 3)

    def keyPressEvent(self, e: QKeyEvent):

        super(ProducedFluidsTable, self).keyPressEvent(e)
        if e.key() == Qt.Key_Delete:
            r = self.currentRow()
            c = self.currentColumn()
            if r == self.rowCount() - 1 and self.rowCount() > 1:
                self.parent().produced_fluids.volumes = np.delete(self.parent().produced_fluids.volumes, (-1), axis=1)
                self.sync_to_data_vols()
            elif r < self.rowCount() - 1 and c == 1:
                self.parent().produced_fluids.volumes[c, r] = np.nan
                self.sync_oil(r)
                self.sync_to_data_vols()
                self.data_changed.emit(r)

        elif e.key() == Qt.Key_Control:
            self.control_pressed = True

        elif e.key() == Qt.Key_C:
            if self.control_pressed:
                indecies = self.selectedIndexes()
                cols = []
                rows = []
                for index in indecies:
                    cols.append(index.column())
                    rows.append(index.row())

                if not all(np.array(cols) == cols[0]):
                    QMessageBox(parent=self.parent(), text='Only copy one column.').exec_()
                    return
                if not cols[0] in [0, 2, 3]:
                    QMessageBox(parent=self.parent(), text='Only copy Aq, Tot, or Oil.').exec_()
                    return

                txt = ''
                for row in rows:
                    txt += self.item(row, cols[0]).text() + '\n'

                QApplication.clipboard().setText(txt)

        elif e.key() == Qt.Key_V:
            if self.control_pressed:
                txt = QApplication.clipboard().text()
                paste_array = self.parse_clipboard(txt)
                if not paste_array.size:
                    QMessageBox(parent=self.parent(), text='Error with paste array.').exec_()
                    return
                indecies = self.selectedIndexes()
                start_row = indecies[0].row()
                col = indecies[0].column()
                if col == 1:
                    QMessageBox(parent=self.parent(), text='Only paste to Aq, Tot, or Oil.').exec_()
                    return
                i = 0
                row_count_changed = False
                self.blockSignals(True)

                for value in list(paste_array.tolist()):
                    row = start_row + i
                    if row > self.rowCount() - 1:
                        self.setRowCount(self.rowCount() + 1)
                        row_count_changed = True
                        # self.rc_changed.emit(self.rowCount())
                    item = self.item(row, col)
                    if item is None:
                        item = QTableWidgetItem('')
                        self.setItem(row, col, item)
                    item.setText(str(np.round(value, 3)))
                    self.blockSignals(False)
                    self.check_item(item)
                    self.blockSignals(True)
                    i += 1
                self.blockSignals(False)

                if self.rowCount() == self.parent().produced_fluids.volumes.shape[1]:
                    row_count_changed = True
                    self.blockSignals(True)
                    self.setRowCount(self.rowCount() + 1)
                    self.blockSignals(False)

                if row_count_changed:
                    self.rc_changed.emit(self.rowCount())

    def parse_clipboard(self, txt: str) -> np.array:

        txt = txt.split('\n')
        if not txt[-1]:
            txt.pop()

        for line in txt:
            if isinstance(line, list):
                return np.array([])

        new_array = []
        for line in txt:

            try:
                value = np.round(float(line), 3)

            except Exception as e:
                QMessageBox(parent=self.parent(), text=str(e)).exec_()
                return np.array([])
            new_array.append(value)

        return np.array(new_array)

    def keyReleaseEvent(self, e: QKeyEvent):
        super(ProducedFluidsTable, self).keyReleaseEvent(e)

        if e.key() == Qt.Key_Control:
            self.control_pressed = False

    def check_item(self, item: QTableWidgetItem):

        c = item.column()
        r = item.row()
        nr = self.rowCount()

        volumes = self.parent().produced_fluids.volumes
        s = volumes.shape

        try:

            value = float(item.text())

            if value < 0.:
                raise ValueError('Value must be non-negative.')
            if c == 1 and r == nr:
                raise ValueError

            value = np.round(value, 3)

        except ValueError:

            print('ValueError detected.')

            if r >= s[1]:
                item = self.takeItem(r, c)
                del item
            else:
                try:
                    if not np.isnan(volumes[c, r]):
                        item.setText(str(volumes[c, r]))
                    else:
                        item.setText('')
                except Exception as e:
                    print(e)

            return

        if r >= s[1]:
            volumes_new = np.resize(volumes, (s[0], s[1] + 1))
            volumes_new[:, :-1] = volumes
            volumes_new[:, -1] = [0., np.nan, 0., 0.]
            volumes = volumes_new

        if c == 1 and (value <= volumes[0, r] or value >= volumes[2, r]):
            value = np.nan

        elif c == 2 and value < volumes[0, r]:
            value = volumes[0, r]

        volumes[c, r] = value

        if c == 2:

            if not np.isnan(volumes[1, r]):
                if value <= volumes[1, r]:
                    volumes[1, r] = np.nan

        elif c == 0:

            if not np.isnan(volumes[1, r]):
                if value >= volumes[1, r]:
                    volumes[1, r] = np.nan

            if value > volumes[2, r]:
                volumes[3, r] = 0.

        self.parent().produced_fluids.volumes = volumes

        if c == 1 or c == 2:
            self.sync_oil(r)
        else:
            self.sync_total(r)

        self.sync_to_data_vols()
        self.data_changed.emit(r)


class SyncTable(QTableWidget):

    def __init__(self, parent: ProducedFluidsView):
        super(SyncTable, self).__init__(parent=parent)

        font = self.font()
        font.setFamily('Courier')
        font.setPointSize(12)
        self.setFont(font)

    def sync_to_data(self, data: np.array):

        self.blockSignals(True)

        if data.shape[1] > 0:

            self.sync_rows()

            for r in range(data.shape[1]):
                for c in range(data.shape[0]):

                    datum = data[c, r]
                    txt = ''
                    if not np.isnan(datum):
                        txt = str(np.round(datum, 3))

                    item = self.item(r, c)
                    if item is None:
                        item = QTableWidgetItem(txt)
                        self.setItem(r, c, item)
                    else:
                        item.setText(txt)

        self.blockSignals(False)

    def sync_rows(self):

        nr = self.parent().table.rowCount()
        self.setRowCount(nr)


class RevisionsTable(SyncTable):

    data_changed = pyqtSignal(int)

    def __init__(self, parent: ProducedFluidsView):
        super(RevisionsTable, self).__init__(parent=parent)

        self.setColumnCount(4)
        self.setHorizontalHeaderLabels(['Aq [mL]', 'M.I. [mL]', 'Tot [mL]', 'Oil [mL]'])

        self.sync_to_data_revisions()

        self.itemChanged.connect(self.check_item)

    def reset_data_at_row(self, r: int):

        pf = self.parent().produced_fluids
        revisions = pf.revisions
        volumes = pf.volumes

        revisions[:, r] = 0.
        if np.isnan(volumes[1, r]):
            revisions[1, r] = np.nan

        self.sync_to_data_revisions()
        self.data_changed.emit(r)

    def sync_to_data_revisions(self):

        self.sync_rows()
        self.sync_to_data(self.parent().produced_fluids.revisions[:, :-1])

    def sync_total(self, r: int):

        pf = self.parent().produced_fluids
        volumes = pf.volumes
        revisions = pf.revisions
        net = volumes + revisions[:, :-1]

        if not np.isnan(net[1, r]):
            revisions[2, r] = np.round(net[1, r] + net[3, r] - volumes[2, r], 3)
        else:
            revisions[2, r] = np.round(net[0, r] + net[3, r] - volumes[2, r], 3)

    def sync_oil(self, r: int):

        pf = self.parent().produced_fluids
        volumes = pf.volumes
        revisions = pf.revisions
        net = volumes + revisions[:, :-1]

        if not np.isnan(net[1, r]):
            revisions[3, r] = np.round(net[2, r] - net[1, r] - volumes[3, r], 3)
        else:
            revisions[3, r] = np.round(net[2, r] - net[0, r] - volumes[3, r], 3)

        if revisions[3, r] == 0.:
            # Correct the potential sign carryover from np.round().
            revisions[3, r] = 0.

    def sync_rows(self):
        super(RevisionsTable, self).sync_rows()

        nr = self.parent().table.rowCount()
        r = self.parent().produced_fluids.revisions
        s = r.shape

        if nr == s[1]:
            return

        if nr > s[1]:
            new_r = np.zeros((s[0], s[1] + 1))
            new_r[:, :s[1]] = r
            new_r[:, -1] = 0.
            new_r[1, -1] = np.nan
            self.parent().produced_fluids.revisions = new_r

        else:
            r = np.delete(r, (-1), axis=1)
            r[:, -1] = 0.
            for col in range(self.columnCount()):
                if self.item(nr - 1, col) is not None:
                    print('deleting item in revisions table.')
                    item = self.takeItem(nr - 1, col)
                    del item

            self.parent().produced_fluids.revisions = r

    def setVisible(self, visible: bool):
        super(RevisionsTable, self).setVisible(visible)
        if visible:
            self.sync_to_data_revisions()

    def keyPressEvent(self, e: QKeyEvent):

        super(RevisionsTable, self).keyPressEvent(e)
        if e.key() == Qt.Key_Delete and self.currentRow() < self.rowCount() - 1:
            self.parent().produced_fluids.revisions[:, self.currentRow()] = 0.
            self.sync_to_data_revisions()

    def check_item(self, item: QTableWidgetItem):

        c = item.column()
        r = item.row()

        if r + 1 == self.rowCount():
            item = self.takeItem(r, c)
            del item
            return

        volumes = self.parent().produced_fluids.volumes

        if c == 1 and np.isnan(volumes[1, r]):
            item.setText('')
            return

        revisions = self.parent().produced_fluids.revisions

        try:
            value = float(item.text())
            value = np.round(value, 3)
            if volumes[c, r] + value < 0.:
                value = -np.round(volumes[c, r], 3)

        except ValueError:
            if not np.isnan(revisions[c, r]):
                item.setText(str(np.round(revisions[c, r], 3)))
            else:
                item.setText('')
            return

        old_value = revisions[c, r]
        revisions[c, r] = value
        net = volumes + revisions[:, :-1]

        if c == 1:

            if net[1, r] <= net[0, r]:
                revisions[1, r] = old_value
            elif net[1, r] >= net[2, r]:
                revisions[1, r] = old_value

        elif c == 0:

            if not np.isnan(net[1, r]):
                if net[0, r] >= net[1, r]:
                    revisions[0, r] = old_value
            elif net[0, r] > net[2, r]:
                revisions[0, r] = old_value

        elif c == 2:

            if not np.isnan(net[1, r]):
                if net[2, r] <= net[1, r]:
                    revisions[2, r] = old_value
            elif net[2, r] < net[0, r]:
                revisions[2, r] = old_value

        elif c == 3:
            pass

        self.parent().produced_fluids.revisions = revisions

        if c == 1 or c == 2 or (c == 0 and np.isnan(volumes[1, r])):
            self.sync_oil(r)
        elif c == 3:
            self.sync_total(r)

        self.sync_to_data_revisions()
        self.data_changed.emit(r)

    def mouseDoubleClickEvent(self, e: QMouseEvent):
        super(RevisionsTable, self).mouseDoubleClickEvent(e)

        indecies = self.selectedIndexes()
        if len(indecies) != 1:
            return

        self.parent().volumes_radios[2].click()
        self.parent().net_table.setCurrentIndex(indecies[0])
        pos = self.verticalScrollBar().value()
        self.parent().net_table.verticalScrollBar().setValue(pos)


class NetTable(SyncTable):

    def __init__(self, parent: ProducedFluidsView):
        super(NetTable, self).__init__(parent=parent)

        self.setRowCount(1)
        self.setColumnCount(4)
        self.setHorizontalHeaderLabels(['Aq [mL]', 'M.I. [mL]', 'Tot [mL]', 'Oil [mL]'])

        self.sync_to_data_net()
        self.itemChanged.connect(self.reset_item)

    def reset_item(self, item: QTableWidgetItem):

        r = item.row()
        c = item.column()

        pf = self.parent().produced_fluids
        datum = pf.volumes[c, r] + pf.revisions[c, r]

        if not np.isnan(datum):
            item.setText(str(np.round(datum, 3)))
        else:
            item.setText('')

    def sync_rows(self):

        nr = self.parent().table.rowCount()
        nr_self = self.rowCount()
        if nr != nr_self:
            super(NetTable, self).sync_rows()
        else:
            return

        if nr < nr_self:
            for col in range(self.columnCount()):
                if self.item(nr - 1, col) is not None:
                    print('deleting item in net table.')
                    item = self.takeItem(nr - 1, col)
                    del item

    def sync_to_data_net(self):

        self.sync_rows()

        pf = self.parent().produced_fluids

        if pf.revisions.shape[1] == pf.volumes.shape[1]:
            data = pf.volumes + pf.revisions
        else:
            data = pf.volumes + pf.revisions[:, :-1]

        self.sync_to_data(data)
        self.set_highlights()

    def set_highlights(self):

        revs = self.parent().produced_fluids.revisions

        for r in range(self.rowCount() - 1):
            for c in range(self.columnCount()):
                color = Qt.black
                if not np.isnan(revs[c, r]):
                    if revs[c, r] != 0.:
                        color = Qt.red
                self.item(r, c).setForeground(color)

    def setVisible(self, visible: bool):
        super(NetTable, self).setVisible(visible)
        if visible:
            self.sync_to_data_net()

    def mouseDoubleClickEvent(self, e: QMouseEvent):
        super(NetTable, self).mouseDoubleClickEvent(e)

        indecies = self.selectedIndexes()
        if len(indecies) != 1:
            return

        self.parent().volumes_radios[1].click()
        self.parent().revisions_table.setCurrentIndex(indecies[0])
        pos = self.verticalScrollBar().value()
        self.parent().revisions_table.verticalScrollBar().setValue(pos)


class METable(SyncTable):

    def __init__(self, parent: ProducedFluidsView):
        super(METable, self).__init__(parent=parent)

        self.setColumnCount(4)
        self.setHorizontalHeaderLabels(['Aq [mL]', 'M.I. [mL]', 'Type', 'Oil Frac.'])

        self.sync_to_data_me()

        self.itemChanged.connect(self.check_item)

    def update_data_at_row(self, r: int):

        pf = self.parent().produced_fluids
        net = pf.volumes[:, r] + pf.revisions[:, r]

        self.item(r, 0).setText(str(np.round(net[0], 3)))

        if np.isnan(net[1]):
            self.item(r, 1).setText('')

            if pf.micro_emulsions[0, r] == 3.:
                pf.micro_emulsions[:, r] = np.nan
                self.item(r, 2).setText('')
                self.item(r, 3).setText('')

        else:
            self.item(r, 1).setText(str(np.round(net[1], 3)))

            if pf.micro_emulsions[0, r] != 3.:
                pf.micro_emulsions[0, r] = 3.
                pf.micro_emulsions[1, r] = 0.5
                self.item(r, 2).setText('III')
                self.item(r, 3).setText('0.5')

    def sync_to_data_me(self):

        self.sync_rows()

        pf = self.parent().produced_fluids
        data = np.zeros(pf.volumes.shape)
        data[:, :] = np.nan
        data[:2, :] = pf.volumes[:2, :] + pf.revisions[:2, :-1]

        if pf.micro_emulsions.shape[1] == pf.volumes.shape[1]:
            data[2:, :] = pf.micro_emulsions[:, :]
        else:
            data[2:, :] = pf.micro_emulsions[:, :-1]

        self.sync_to_data(data)

    def sync_to_data(self, data: np.array):
        super(METable, self).sync_to_data(data)

        for row in range(self.rowCount() - 1):
            item = self.item(row, 2)
            if item.text() in ['1', '1.0']:
                item.setText('I')
            elif item.text() in ['2', '2.0']:
                item.setText('II')
            elif item.text() in ['3', '3.0']:
                item.setText('III')

    def setVisible(self, visible: bool):
        super(METable, self).setVisible(visible)
        if visible:
            self.sync_to_data_me()

    def sync_rows(self):
        super(METable, self).sync_rows()

        nr = self.parent().table.rowCount()
        me = self.parent().produced_fluids.micro_emulsions
        s = me.shape

        if nr == s[1]:
            return

        if nr > s[1]:
            new_me = np.zeros((s[0], s[1] + 1))
            new_me[:, :s[1]] = me
            new_me[:, -1] = np.nan
            self.parent().produced_fluids.micro_emulsions = new_me
            self.sync_to_data_me()

        else:
            me = np.delete(me, (-1), axis=1)
            me[:, -1] = np.nan
            for col in range(self.columnCount()):
                if self.item(nr - 1, col) is not None:
                    item = self.takeItem(nr - 1, col)
                    del item

            self.parent().produced_fluids.micro_emulsions = me

    def check_item(self, item: QTableWidgetItem):

        c = item.column()
        r = item.row()

        if r + 1 == self.rowCount():
            item = self.takeItem(r, c)
            del item
            return

        volumes = self.parent().produced_fluids.volumes
        revisions = self.parent().produced_fluids.revisions
        net = volumes + revisions[:, :-1]
        me = self.parent().produced_fluids.micro_emulsions

        if c == 0 or c == 1:
            txt = ''
            if not (c == 1 and np.isnan(net[c, r])):
                txt = str(np.round(net[c, r], 3))
            item.setText(txt)
            return
        elif c == 2:
            if not np.isnan(net[1, r]):
                item.setText('III')
                return
        elif c == 3 and np.isnan(me[0, r]):
            item.setText('')
            return

        txt = item.text()
        me_dict = {'I': 1., 'II': 2., 'III': 3.}

        try:
            if c == 2:
                if txt.upper() in me_dict.keys():
                    value = me_dict[txt.upper()]
                else:
                    value = float(txt)
                    if value not in me_dict.values():
                        raise ValueError

                if np.isnan(volumes[1, r]) and value == 3.:
                    raise ValueError

                old_value = me[0, r]
                if old_value != value or np.isnan(old_value):
                    if value == 1:
                        me[1, r] = 0.
                    elif value == 2:
                        me[1, r] = 1.
                    else:
                        me[1, r] = 0.5

                    self.blockSignals(True)
                    self.item(r, 3).setText(str(np.round(me[1, r], 3)))
                    self.blockSignals(False)

                new_txt = 'I'
                if value == 2.:
                    new_txt = 'II'
                elif value == 3.:
                    new_txt = 'III'

            else:
                value = np.round(float(txt), 3)
                if not 0. <= value <= 1.:
                    raise ValueError

                new_txt = str(value)

            item.setText(new_txt)

        except ValueError:

            if not np.isnan(me[c - 2, r]):
                if c == 2:
                    txt = 'I'
                    if me[0, r] == 2.:
                        txt = 'II'
                    elif me[0, r] == 3:
                        txt = 'III'
                else:
                    txt = str(np.round(me[c - 2, r], 3))

                item.setText(txt)
            else:
                item.setText('')
            return

        me[c - 2, r] = value

    def keyPressEvent(self, e: QKeyEvent):

        if e.key() == Qt.Key_Delete:

            r = self.currentRow()
            c = self.currentColumn()
            me = self.parent().produced_fluids.micro_emulsions
            volumes = self.parent().produced_fluids.volumes

            if c == 2:

                if not np.isnan(volumes[1, r]):
                    pass
                else:
                    me[:, r] = np.nan
                    self.sync_to_data_me()

            elif c == 3:

                if not np.isnan(me[0, r]):

                    if me[0, r] == 1.:
                        me[1, r] = 0.
                    elif me[0, r] == 2.:
                        me[1, r] = 1.
                    else:
                        me[1, r] = 0.5

                else:
                    me[1, r] = np.nan

                self.sync_to_data_me()

        super(METable, self).keyPressEvent(e)


class FinalValuesTable(SyncTable):

    def __init__(self, parent: ProducedFluidsView):
        super(FinalValuesTable, self).__init__(parent=parent)

        self.setColumnCount(3)
        self.setHorizontalHeaderLabels(['Aq [mL]', 'Tot [mL]', 'Oil [mL]'])

        self.sync_to_data_final()

        # self.itemChanged.connect(self.check_item)

    def sync_to_data_final(self):
        print('sync_to_data_final')

        self.sync_rows()

        pf = self.parent().produced_fluids
        print('getting final oil and water')
        final = pf.final_oil_and_water()
        print('got final oil and water')

        self.sync_to_data(final)

    def setVisible(self, visible: bool):
        super(FinalValuesTable, self).setVisible(visible)
        if visible:
            self.sync_to_data_final()


class EffluentMeasurementsTable(QTableWidget):

    def __init__(self, parent: ProducedFluidsView):
        super(EffluentMeasurementsTable, self).__init__(parent=parent)

        sf = parent.width() / 1010.
        self.setFixedWidth(int(sf * 550))
        font = self.font()
        font.setFamily('Courier')
        font.setPointSize(12)
        self.setFont(font)

        self.setRowCount(1)
        measurements = parent.produced_fluids.measurements
        s = measurements.shape
        self.setColumnCount(s[0])
        self.setHorizontalHeaderLabels(['pH', 'R.I.', 'ORP [mV]', chr(951) + u' [cP]', 'A520'])

        if s[1] > 0:

            self.setRowCount(s[1] + 1)

            for r in range(s[1]):
                for c in range(s[0]):
                    if not np.isnan(measurements[c, r]):
                        self.setItem(r, c, QTableWidgetItem(str(measurements[c, r])))
                    else:
                        self.setItem(r, c, QTableWidgetItem(''))

        self.itemChanged.connect(self.check_item)

    def sync_rows(self, nr: int):

        nr = self.parent().table.rowCount()
        self.setRowCount(nr)
        m = self.parent().produced_fluids.measurements
        s = m.shape

        if nr == s[1]:
            return

        if nr > s[1]:
            new_m = np.zeros((s[0], s[1] + 1))
            new_m[:, :s[1]] = m
            new_m[:, -1] = np.nan
            self.parent().produced_fluids.measurements = new_m
            m = new_m
        else:
            m = np.delete(m, (-1), axis=1)
            m[:, -1] = np.nan
            for col in range(self.columnCount()):
                if self.item(nr - 1, col) is not None:
                    item = self.takeItem(nr - 1, col)
                    del item

            self.parent().produced_fluids.measurements = m

    def check_item(self, item: QTableWidgetItem):

        r = item.row()
        c = item.column()
        m = self.parent().produced_fluids.measurements
        s = m.shape

        if r + 1 == self.rowCount():
            item = self.takeItem(r, c)
            del item
            return

        try:
            value = float(item.text())
        except ValueError:
            if not np.isnan(self.parent().produced_fluids.measurements[c, r]):
                item.setText(str(self.parent().produced_fluids.measurements[c, r]))
            else:
                item.setText('')
            return

        if r >= s[1]:
            new_m = np.resize(m, (s[0], s[1] + 1))
            new_m[:, :-1] = m
            new_m[:, -1] = np.nan
            m = new_m

        m[c, r] = value
        self.parent().produced_fluids.measurements = m

    def keyPressEvent(self, e: QKeyEvent):

        if e.key() == Qt.Key_Delete:

            r = self.currentRow()
            c = self.currentColumn()

            self.takeItem(r, c)

            self.parent().produced_fluids.measurements[c, r] = np.nan

        super(EffluentMeasurementsTable, self).keyPressEvent(e)


class DilutionMeasurementsTable(SyncTable):

    def __init__(self, parent: ProducedFluidsView):
        super(DilutionMeasurementsTable, self).__init__(parent=parent)

        sf = parent.width() / 1010.
        self.setFixedWidth(int(sf * 550))

    def check_item(self, item: QTableWidgetItem):

        r = item.row()
        c = item.column()

        if r + 1 == self.rowCount():
            item = self.takeItem(r, c)
            del item
            return


class IonsTable(DilutionMeasurementsTable):

    def __init__(self, parent: ProducedFluidsView):
        super(IonsTable, self).__init__(parent=parent)

        ions = parent.produced_fluids.ions
        s = ions.shape
        self.setColumnCount(s[0])
        self.sync_rows()
        self.setHorizontalHeaderLabels(['Samp [g]', 'Dil [g]', 'Tot [g]', 'Na+ [ppm]', 'K+ [ppm]', 'Mg++ [ppm]',
                                        'Ca++ [ppm]', 'Cl- [ppm]', 'Br- [ppm]', 'SO4-- [ppm]'])
        for i in range(3, 10):
            self.setColumnWidth(i, 120)

        self.itemChanged.connect(self.check_item)
        self.horizontalHeader().sectionClicked.connect(self.set_tracer_target)
        self.sync_to_data_ions()

    def sync_to_data_ions(self):

        data = self.parent().produced_fluids.ions
        self.sync_to_data(data)

    def check_item(self, item: QTableWidgetItem):
        print('checking item.')
        super(IonsTable, self).check_item(item)

        r = item.row()
        c = item.column()
        print(r, c)
        ions = self.parent().produced_fluids.ions
        # s = ions.shape

        try:
            value = float(item.text())
            if value < 0.:
                raise ValueError

        except ValueError:
            if not np.isnan(ions[c, r]):
                item.setText(str(np.round(ions[c, r], 3)))
            else:
                item.setText('')
            return

        print('setting value = {}', value)
        print(ions)
        ions[c, r] = value
        print('item set.')

        if c == 0 or c == 1:

            ions[2, r] = ions[0, r] + ions[1, r]
            if not np.isnan(ions[2, r]):
                ions[2, r] = np.round(ions[2, r], 3)

        elif c == 2:

            ions[1, r] = ions[2, r] - ions[0, r]
            if not np.isnan(ions[1, r]):
                ions[1, r] = np.round(ions[1, r], 3)

        print('syncing data to ions.')
        self.sync_to_data_ions()

    def keyPressEvent(self, e: QKeyEvent):
        super(IonsTable, self).keyPressEvent(e)

        if e.key() == Qt.Key_Delete:

            ions = self.parent().produced_fluids.ions
            c = self.currentColumn()
            r = self.currentRow()

            if c == 0:
                ions[:3, r] = np.nan
                self.sync_to_data_ions()

            elif c == 1 or c == 2:
                pass

            elif c > 2:
                ions[c, r] = np.nan
                self.sync_to_data_ions()

    def set_tracer_target(self, *args):

        self.parent().tracer_target = args[0]


class AqueousTracer:

    def __init__(self, name: str, cum_vols: np.array, cum_aq_vols: np.array, aq_conc: np.array):

        self.name = name
        self.cv = cum_vols
        self.cv_aq = cum_aq_vols
        self.tracer = aq_conc
        self._ipv = 0.

        non_nan = np.where(~np.isnan(self.tracer))

        if not non_nan[0].tolist():
            return

        self.cv_original = self.cv.copy()
        self.cv_aq_original = self.cv_aq.copy()
        s = self.cv_original.shape
        cv2 = np.zeros((s[0] + 1,))
        cv2[1:] = self.cv_original
        cv2w = np.zeros((s[0] + 1,))
        cv2w[1:] = self.cv_aq_original
        self.volumes_original = np.diff(cv2)
        self.aq_volumes_original = np.diff(cv2w)
        self.fw_original = np.divide(self.aq_volumes_original, self.volumes_original)

        self.cv = self.cv[non_nan]
        self.cv_aq = self.cv_aq[non_nan]

        self.tracer = self.tracer[non_nan]
        self.non_nan = non_nan
        self.n_min = np.nan
        self.n_max = np.nan

        s = self.cv.shape
        cv2 = np.zeros((s[0] + 1,))
        cv2[1:] = self.cv
        cv2w = np.zeros((s[0] + 1,))
        cv2w[1:] = self.cv_aq
        self.volumes = np.diff(cv2)
        self.aq_volumes = np.diff(cv2w)
        self.fw = np.divide(self.aq_volumes, self.volumes)

        self.excluded = []
        for i in range(self.tracer.size):
            self.excluded.append(False)

    def get_ipv(self) -> float:

        return self._ipv

    def set_ipv(self, val: float):

        if val < 0.:
            raise ValueError

        self._ipv = val

    def exclude_include_point(self, index: int):

        if self.excluded[index]:
            self.excluded[index] = False
        else:
            self.excluded[index] = True

    def is_point_excluded(self, index: int):

        return self.excluded[index]

    def set_normalization_min(self, nm: float):

        self.n_min = nm

    def set_normalization_max(self, nm: float):

        self.n_max = nm

    def adjust_min_max(self, n_min: float, n_max: float):

        if not np.isnan(self.n_min):
            n_min = self.n_min

        if not np.isnan(self.n_max):
            n_max = self.n_max

        return n_min, n_max

    def get_adjusted_tracer(self) -> np.array:

        return np.multiply(self.fw_original[self.non_nan], self.tracer)

    def get_normalized_tracer(self) -> np.array:

        n_min = self.tracer.min()
        n_max = self.tracer.max()
        n_min, n_max = self.adjust_min_max(n_min, n_max)

        return (self.tracer - n_min) / (n_max - n_min)

    def get_inc_normalized_tracer(self) -> np.array:

        old_tracer = self.tracer.copy()
        new_tracer = old_tracer[np.where(np.array(self.excluded) == 0)]

        n_min = new_tracer.min()
        n_max = new_tracer.max()
        n_min, n_max = self.adjust_min_max(n_min, n_max)

        return (new_tracer - n_min) / (n_max - n_min)

    def get_normalized_adjusted_tracer(self) -> np.array:

        return np.multiply(self.fw_original[self.non_nan], self.get_normalized_tracer())

    def get_inc_normalized_adj_tracer(self) -> np.array:

        fw = self.fw_original.copy()[self.non_nan]
        return np.multiply(fw[np.where(np.array(self.excluded) == 0)], self.get_inc_normalized_tracer())

    def get_inc_cv(self) -> np.array:

        return self.cv[np.where(np.array(self.excluded) == 0)]

    def get_inc_cv_aq(self) -> np.array:

        return self.cv_aq[np.where(np.array(self.excluded) == 0)]

    @staticmethod
    def assemble_staircase(cv: list, y_data: list, final_x_val: float):

        s_x_data = [0.]
        s_y_data = [y_data[0]]

        for i in range(len(cv)):
            if i == 0:
                s_x_data.append(0.)
                s_y_data.append(y_data[i])
            else:
                s_x_data.append(cv[i - 1])
                s_y_data.append(y_data[i])

            s_x_data.append(cv[i])
            s_y_data.append(y_data[i])

        s_x_data.append(final_x_val)
        s_y_data.append(s_y_data[-1])

        return s_x_data, s_y_data

    def get_original_staircase(self):

        cv = self.cv_original.copy().tolist()
        cv_aq = self.cv_aq.copy().tolist()
        y_data = self.get_normalized_tracer().tolist()

        f = interp1d(cv_aq, y_data, bounds_error=False, fill_value=np.nan)
        y_data_interp = f(self.cv_aq_original)
        y_data_interp = np.multiply(y_data_interp, self.fw_original)

        s_x_data, s_y_data = AqueousTracer.assemble_staircase(cv, y_data_interp, self.cv_original[-1])

        return s_x_data, s_y_data

    def get_staircase(self):

        cv = self.cv_original.copy().tolist()
        cv_aq = self.cv_aq_original.copy().tolist()
        cv_aq_inc = self.get_inc_cv_aq().tolist()
        y_data_inc = self.get_inc_normalized_tracer().tolist()

        f = interp1d(cv_aq_inc, y_data_inc, bounds_error=False, fill_value=np.nan)
        y_data_interp = f(np.array(cv_aq))
        y_data_interp = np.multiply(y_data_interp, self.fw_original)

        s_x_data, s_y_data = AqueousTracer.assemble_staircase(cv, y_data_interp.tolist(), self.cv_original[-1])

        return s_x_data, s_y_data

    def get_original_center_points(self):

        x_data = self.cv_original.copy() - 0.5 * self.volumes_original.copy()
        x_data = x_data[self.non_nan]
        y_data = self.get_normalized_adjusted_tracer()

        return x_data.tolist(), y_data.tolist()

    def get_center_points(self):

        x_data = self.cv_original.copy() - 0.5 * self.volumes_original.copy()
        x_data = x_data[self.non_nan]
        cv_aq = self.cv_aq.copy()
        cv_aq_inc = self.get_inc_cv_aq()
        y_data_inc = self.get_inc_normalized_tracer()

        f = interp1d(cv_aq_inc, y_data_inc, bounds_error=False, fill_value=np.nan)
        y_data_interp = f(cv_aq)
        y_data_interp = np.multiply(y_data_interp, self.fw_original[self.non_nan])

        return x_data.tolist(), y_data_interp.tolist()

    def get_all_center_points(self):

        x_data = self.cv_original.copy() - 0.5 * self.volumes_original.copy()
        cv_aq = self.cv_aq_original.copy()
        cv_aq_inc = self.get_inc_cv_aq()
        y_data_inc = self.get_inc_normalized_tracer()

        f = interp1d(cv_aq_inc, y_data_inc, bounds_error=False, fill_value=np.nan)
        y_data_interp = f(cv_aq)
        y_data_interp = np.multiply(y_data_interp, self.fw_original)

        return x_data.tolist(), y_data_interp.tolist()

    def get_interpolated_tracer(self) -> np.array:

        cv_aq = self.cv_aq_original.copy()
        cv_aq_inc = self.get_inc_cv_aq()
        y_data = self.tracer.copy()
        y_data_inc = []

        for yd, ex in zip(y_data, self.excluded):
            if not ex:
                y_data_inc.append(yd)

        y_data_inc = np.array(y_data_inc)

        f = interp1d(cv_aq_inc, y_data_inc, bounds_error=False, fill_value=np.nan)
        y_data_interp = f(cv_aq)

        n_min = np.min(y_data_inc)
        n_max = np.max(y_data_inc)

        n_min, n_max = self.adjust_min_max(n_min, n_max)
        no_nan = np.where(~np.isnan(y_data_interp))
        if no_nan[0][0] > 0:
            y_data_interp[:no_nan[0][0]] = n_min
        if no_nan[0][-1] < y_data_interp.size - 1:
            y_data_interp[no_nan[0][-1] + 1:] = n_max

        return y_data_interp

    def is_tracer_area_above(self) -> bool:

        _, cp = self.get_all_center_points()

        if cp[-1] > cp[0]:
            return True

        return False

    def get_area_above_below_tracer(self, above: bool) -> float:

        _, cp = self.get_all_center_points()
        cp = np.array(cp)
        non_nans = np.where(~np.isnan(cp))
        cp = cp[non_nans]
        v_non_nan = self.volumes_original.copy()
        v_non_nan = v_non_nan[non_nans]

        if above:
            return np.sum(np.multiply(1. - cp, v_non_nan)) + 0.
        else:
            return np.sum(np.multiply(cp, v_non_nan)) + 0.

    def get_aq_pv_apparent(self) -> float:

        _, cp = self.get_all_center_points()
        cp = np.array(cp)
        non_nans = np.where(~np.isnan(cp))
        cv_non_nan = self.cv_original.copy()
        cv_non_nan = cv_non_nan[non_nans]
        cv_aq_non_nan = self.cv_aq_original.copy()
        cv_aq_non_nan = cv_aq_non_nan[non_nans]
        v_non_nan = self.volumes_original.copy()
        v_non_nan = v_non_nan[non_nans]
        additional_oil = self.cv_original[-1] - self.cv_aq_original[-1] - cv_non_nan[-1] + cv_aq_non_nan[-1]
        previous_fluid = cv_non_nan[0] - v_non_nan[0]
        all_oil = self.cv_original[-1] - self.cv_aq_original[-1]

        above = self.is_tracer_area_above()
        tracer_area = self.get_area_above_below_tracer(above)

        if above:
            return tracer_area + additional_oil + previous_fluid
        else:
            return tracer_area + all_oil

    def get_aq_pv(self, true_aq_pv: float) -> float:

        return self.get_aq_pv_apparent() + true_aq_pv * self._ipv

    def get_retention(self, true_aq_pv: float) -> float:

        return self.get_aq_pv(true_aq_pv) - true_aq_pv

    def get_d_factor(self, true_aq_pv: float) -> float:

        return self.get_retention(true_aq_pv) / true_aq_pv


class AqueousTracerView(QWidget):

    staircase_updated = pyqtSignal(list)
    link_tracer = pyqtSignal(bool)

    def __init__(self, parent, tracer: AqueousTracer):
        super(AqueousTracerView, self).__init__(parent=parent)

        self.tracer_object = None

        self.setLayout(QVBoxLayout())

        self.link_button = QPushButton(parent=self, text='Link')
        self.link_button.clicked.connect(self.link_button_clicked)
        font = self.link_button.font()
        font.setPointSize(12)
        self.link_button.setFont(font)

        self.layout().addWidget(self.link_button)

        self.chart = QChart(flags=Qt.WindowFlags())
        self.chart.legend().hide()
        self.view = flood.CoreFloodChartView(self.chart, self)
        self.x_axis = QValueAxis()
        self.y_axis = QValueAxis()
        self.chart.setAxisX(self.x_axis)
        self.chart.setAxisY(self.y_axis)
        font = QFont()
        font.setPointSize(11)
        font.setFamily('Courier')
        title_font = QFont()
        title_font.setPointSize(12)
        title_font.setBold(True)
        self.x_axis.setLabelsFont(font)
        self.y_axis.setLabelsFont(font)
        self.x_axis.setTitleFont(title_font)
        self.y_axis.setTitleFont(title_font)
        self.x_axis.setTitleText('Volume Injected [mL]')

        self.layout().addWidget(self.view)

        self.data = []
        self.data_color = Qt.black
        self.original_staircase = QLineSeries()
        self.staircase = QLineSeries()
        self.x_max = 1.
        self.add_series_to_plot(self.original_staircase)
        self.original_staircase.setColor(Qt.gray)
        self.add_series_to_plot(self.staircase)

        self.y_axis.setMax(1.)

        self.set_tracer(tracer)
        self.set_data_color(Qt.black)
        self.interactive = True
        self.show()

    def link_button_clicked(self):

        if not self.interactive:
            return

        if self.link_button.text() == 'Link':
            self.link_tracer.emit(True)
            self.link_button.setText('Un-link')
        else:
            self.link_tracer.emit(False)
            self.link_button.setText('Link')

    def show_data(self, visible: bool):

        self.original_staircase.setVisible(visible)
        self.staircase.setVisible(visible)
        for datum in self.data:
            datum.setVisible(visible)
        self.interactive = visible

    def set_data_color(self, color):

        self.data_color = color
        self.staircase.setColor(color)
        for datum, ex in zip(self.data, self.tracer_object.excluded):
            if not ex:
                datum.setColor(color)

    def add_series_to_plot(self, series):

        self.chart.addSeries(series)
        series.attachAxis(self.x_axis)
        series.attachAxis(self.y_axis)
        series.setUseOpenGL(False)

    def plot_data(self):

        x_data, y_data = self.tracer_object.get_original_center_points()
        self.x_max = self.tracer_object.cv_original.copy()[-1]

        for xd, yd, exclude in zip(x_data, y_data, self.tracer_object.excluded):
            data_line = QScatterSeries()
            if not exclude:
                data_line.setColor(Qt.black)
            else:
                data_line.setColor(Qt.gray)
            data_line.setMarkerShape(QScatterSeries.MarkerShapeCircle)
            data_line.setMarkerSize(10.)
            point = utils.series_to_polyline([xd], [yd])
            self.add_series_to_plot(data_line)
            data_line.append(point)
            self.data.append(data_line)

        s_x_data, s_y_data = self.tracer_object.get_original_staircase()

        stair = utils.series_to_polyline(s_x_data, s_y_data)
        self.original_staircase.clear()
        self.original_staircase.append(stair)
        self.staircase.clear()
        self.staircase.append(stair)

    def update_staircase(self):

        s_x_data, s_y_data = self.tracer_object.get_staircase()

        stair = utils.series_to_polyline(s_x_data, s_y_data)
        self.staircase.clear()
        self.staircase.append(stair)

        xds, yds = self.tracer_object.get_center_points()
        _, ydos = self.tracer_object.get_original_center_points()

        for datum, xd, yd, ydo in zip(self.data, xds, yds, ydos):
            if datum.color() == Qt.black:
                point = utils.series_to_polyline([xd], [yd])
            else:
                point = utils.series_to_polyline([xd], [ydo])
            datum.replace(point)

        self.staircase_updated.emit(self.tracer_object.excluded)

    def clear_plot(self):

        self.original_staircase.clear()
        self.staircase.clear()
        for i in range(len(self.data)):
            datum = self.data.pop()
            datum.clear()
            del datum

    def set_tracer(self, tracer: AqueousTracer):

        self.tracer_object = tracer
        self.clear_plot()
        self.plot_data()
        self.update_staircase()
        self.y_axis.setTitleText('Normalized ' + tracer.name)
        self.x_axis.setMin(0.)
        self.x_axis.setMax(self.x_max)

    def mousePressEvent(self, a0: QMouseEvent):

        if a0.button() != 1 or not self.interactive:
            return

        parent_point = QPoint(a0.x(), a0.y())
        point = self.view.mapFromParent(parent_point)
        scene_point = self.view.mapToScene(point)
        chart_point = self.chart.mapFromScene(scene_point)
        value_point = self.chart.mapToValue(chart_point)

        if value_point.x() > self.x_axis.max() or value_point.x() < self.x_axis.min():
            return

        if value_point.y() > self.y_axis.max() or value_point.y() < self.y_axis.min():
            return

        p1 = self.chart.mapToValue(QPointF(0., 0.))
        p2 = self.chart.mapToValue(QPointF(100., 0.))
        p3 = self.chart.mapToValue(QPointF(0., 100.))
        mx = (p2.x() - p1.x()) / 100.
        my = (p3.y() - p1.y()) / 100.

        x_data, y_data = self.tracer_object.get_center_points()
        _, y_data_o = self.tracer_object.get_original_center_points()
        i = 0
        for ex, ydo in zip(self.tracer_object.excluded, y_data_o):
            if ex:
                y_data[i] = ydo
            i += 1
        x_data = np.array(x_data)
        y_data = np.array(y_data)
        x_data_mapped = (x_data - p1.x()) / mx
        y_data_mapped = (y_data - p1.y()) / my

        index = utils.find_nearest(x_data_mapped, y_data_mapped, chart_point.x(), chart_point.y())

        self.tracer_object.exclude_include_point(index)
        if self.tracer_object.excluded[index]:
            self.data[index].setColor(Qt.gray)
        else:
            self.data[index].setColor(self.data_color)

        self.update_staircase()


class PorosityTool(QDialog):

    def __init__(self, parent: ProducedFluidsView):
        super(PorosityTool, self).__init__(parent=parent)

        self.setWindowTitle('Porosity Tool')
        layout = QVBoxLayout()
        h_layout = QHBoxLayout()
        h_layout_2 = QHBoxLayout()
        self.setLayout(layout)
        self.porosity_text = QLabel(parent=self, text='Aqueous Pore Volume [mL]: ')
        font = self.porosity_text.font()
        font.setPointSize(12)
        self.porosity_text.setFont(font)
        h_layout.addWidget(self.porosity_text)
        self.send_porosity_button = QPushButton(parent=self, text='Set as total porosity')
        self.send_porosity_button.setFont(font)
        self.send_porosity_button.clicked.connect(self.set_experiment_porosity)

        h_layout_2.addWidget(self.send_porosity_button)
        try:
            if parent.parent().flood.use_for_ref:
                txt = 'Do not use for ref'
            else:
                txt = 'Use for ref'
        except Exception as e:
            print(e)
            txt = 'Use for ref'
            parent.parent().flood.__dict__['use_for_ref'] = False
        self.use_for_ref_button = QPushButton(parent=self, text=txt)
        self.use_for_ref_button.setFont(font)
        self.use_for_ref_button.clicked.connect(self.use_for_ref)

        self.dead_volume_label = QLabel(parent=self, text='- Dead Volume [mL]: ')
        self.dead_volume_label.setFont(font)
        self.dead_volume_edit = QLineEdit(parent=self)
        self.dead_volume_edit.setFont(font)
        self.sync_dead_volume_edit()
        self.net_volume_text = QLabel(parent=self, text='= mL')
        self.net_volume_text.setFont(font)
        h_layout.addWidget(self.dead_volume_label)
        h_layout.addWidget(self.dead_volume_edit)
        h_layout.addWidget(self.net_volume_text)
        h_layout_2.addWidget(self.use_for_ref_button)
        layout.addLayout(h_layout)
        layout.addLayout(h_layout_2)

        self.dead_volume_edit.editingFinished.connect(self.check_dead_volume_edit)

        sf = parent.width() / 1010.
        self.setFixedSize(int(sf * 550), int(sf * 550))

        if parent.parent().flood.tracer_object is None:
            if parent.tracer_table == parent.measurements_table:
                self.tracer = parent.produced_fluids.measurements[1, :]
            else:
                sample = parent.produced_fluids.ions[0, :]
                total = parent.produced_fluids.ions[2, :]
                tracer = parent.produced_fluids.ions[parent.tracer_target, :]
                self.tracer = np.multiply(np.divide(total, sample), tracer)

            non_nan = np.where(~np.isnan(self.tracer))

            if not non_nan[0].tolist():
                QMessageBox(parent=parent, text='No non-nan data available.').exec_()
                return

            if parent.tracer_table == parent.measurements_table:
                name = 'R.I.'
            else:
                if parent.tracer_target < 3:
                    return
                txt = parent.ions_table.horizontalHeaderItem(parent.tracer_target).text()
                txt = txt.split('[')[0].strip()
                name = txt
            self.cv = parent.produced_fluids.cumulative_volume()
            self.cv_aq = parent.produced_fluids.cumulative_aq_volume()
            self.tracer_object = AqueousTracer(name, self.cv.copy(), self.cv_aq.copy(), self.tracer.copy())
            self.send_porosity_button.setEnabled(False)
            self.use_for_ref_button.setEnabled(False)
        else:
            self.tracer_object = parent.parent().flood.tracer_object

        self.tracer_view = AqueousTracerView(parent=self, tracer=self.tracer_object)
        if parent.parent().flood.tracer_object is not None:
            self.tracer_view.link_button.setText('Un-link')
        self.tracer_view.staircase_updated.connect(self.update_pv_and_nv)
        self.tracer_view.link_tracer.connect(self.link_button_clicked)
        layout.addWidget(self.tracer_view)

        self.pv = np.nan
        self.update_pv_and_nv([])
        self.show()

    def update_porosity_text(self):

        self.porosity_text.setText('Aqueous Pore Volume [mL]: {:.1f}'.format(self.pv))

    def check_dead_volume_edit(self):

        txt = self.dead_volume_edit.text()

        try:
            value = float(txt)
            if value < 0.:
                raise ValueError
            value = np.round(value, 1)

        except ValueError:
            value = self.parent().produced_fluids.dead_volume

        self.parent().produced_fluids.dead_volume = value
        self.sync_dead_volume_edit()
        self.update_net_volume()

    def sync_dead_volume_edit(self):

        self.dead_volume_edit.setText(str(np.round(self.parent().produced_fluids.dead_volume, 1)))

    def update_net_volume(self):

        dv = self.parent().produced_fluids.dead_volume
        self.net_volume_text.setText('= {:.1f} mL'.format(self.pv - dv))

    def calc_porosity(self):

        return self.tracer_object.get_aq_pv_apparent()

    def link_button_clicked(self, link: bool):

        if link:
            self.update_flood_pv()
            self.set_flood_tracer_object()
            self.send_porosity_button.setEnabled(True)
            self.use_for_ref_button.setEnabled(True)
        else:
            self.clear_flood_tracer_object()
            self.send_porosity_button.setEnabled(False)
            if self.parent().parent().flood.use_for_ref:
                self.use_for_ref()
            self.use_for_ref_button.setEnabled(False)

    def update_flood_pv(self):

        nv = self.pv - self.parent().produced_fluids.dead_volume
        self.parent().parent().flood.pore_volume = nv
        self.parent().parent().flood.tracer_object = self.tracer_object

    def set_flood_tracer_object(self):

        self.parent().parent().flood.tracer_object = self.tracer_object

    def clear_flood_tracer_object(self):

        self.parent().parent().flood.tracer_object = None

    def set_experiment_porosity(self):

        pore_volume = self.pv
        dead_volume = self.parent().produced_fluids.dead_volume
        net_volume = pore_volume - dead_volume
        experiment = self.parent().parent().linked_widget.experiment
        core = experiment.core
        bulk_volume = core.my_bulk_volume()
        porosity = net_volume / bulk_volume
        experiment.petro_parameters['phi'][1] = porosity
        experiment.petro_parameters['phi'][2] = self.parent().parent().flood

    def use_for_ref(self):

        fl = self.parent().parent().flood
        if fl.use_for_ref:
            fl.use_for_ref = False
            self.use_for_ref_button.setText('Use for ref')
        else:
            fl.use_for_ref = True
            self.use_for_ref_button.setText('Do not use for ref')

    def update_pv_and_nv(self, _):

        self.pv = self.calc_porosity()
        self.update_porosity_text()
        self.update_net_volume()

    def closeEvent(self, _):

        try:
            self.parent().p_tool = None
        except Exception as e:
            print(e)


class RetentionTool(QDialog):

    def __init__(self, parent: ProducedFluidsView):
        super(RetentionTool, self).__init__(parent=parent)

        sf = parent.width() / 1010.
        self.setFixedSize(int(sf * 900), int(sf * 500))
        self.setWindowTitle('Retention Tool')

        lyt = QHBoxLayout()
        self.setLayout(lyt)

        cg_widget = QWidget(self)
        cg_lyt = QGridLayout()
        cg_widget.setLayout(cg_lyt)
        self.va_combo = QComboBox(self)
        self.va_combo.addItems(['Viscosity', 'Absorbance'])
        self.va_combo.currentIndexChanged.connect(self.change_palette)
        cg_lyt.addWidget(self.va_combo, 0, 1, 1, 1)

        lyt.addWidget(cg_widget)
        self.cg_widget = cg_widget

        fl = parent.parent().flood

        self.viscosity_tracer_object = fl.effluent.cp_tracer_object
        self.abs_tracer_object = fl.effluent.abs_tracer_object

        x = parent.produced_fluids.cumulative_volume()
        x_aq = parent.produced_fluids.cumulative_aq_volume()
        y = parent.produced_fluids.measurements[3, :]
        y2 = parent.produced_fluids.measurements[4, :]
        non_nans = np.where(~np.isnan(y))
        non_nans2 = np.where(~np.isnan(y2))

        v_data_exists = True
        if non_nans[0].size == 0:
            v_data_exists = False

        a_data_exists = True
        if non_nans2[0].size == 0:
            a_data_exists = False

        if isinstance(fl, flood.MultiCoreFlood) and (not a_data_exists and self.abs_tracer_object is None):
            raise('Error: Mult-Flood requires absorbance data')

        if not v_data_exists and not a_data_exists and self.viscosity_tracer_object is None and \
                self.abs_tracer_object is None:
            raise('Error: must provide viscosity or absorbance data')

        bfs = fl.experiment.binned_floods(fl)
        i = bfs.index(fl)
        ppm_prev = 0.
        viscosity_prev = bfs[i].get_last_fluid().specific_fluid.brine.get_viscosity(fl.temperature)
        abs_prev = np.nan

        if i > 0:
            fluid_prev = bfs[i - 1].get_last_fluid()
            if fluid_prev.get_fluid_type() == FluidType.POLYMER_SOLUTION:
                ppm_prev = fluid_prev.specific_fluid.polymer_solution.concentration
                viscosity_prev = fluid_prev.get_viscosity(fl.temperature, parent.produced_fluids.gamma)
            elif fluid_prev.get_fluid_type() == FluidType.BRINE:
                viscosity_prev = fluid_prev.brine.get_viscosity(fl.temperature)

        self.ppm_prev = ppm_prev
        if not np.isnan(parent.produced_fluids.retention):
            self.ppm_prev = parent.produced_fluids.retention_ppm_min

        self.viscosity_prev = viscosity_prev
        self.abs_prev = abs_prev

        self.tracer_view = None

        if v_data_exists and self.viscosity_tracer_object is None:
            self.viscosity_tracer_object = AqueousTracer('cP', x, x_aq, y)
            self.viscosity_tracer_object.set_normalization_min(viscosity_prev)
            self.tracer_view = AqueousTracerView(self, self.viscosity_tracer_object)
            self.tracer_view.set_data_color(Qt.blue)
        elif self.viscosity_tracer_object is not None:
            self.tracer_view = AqueousTracerView(self, self.viscosity_tracer_object)
            self.tracer_view.set_data_color(Qt.blue)

        if a_data_exists and self.abs_tracer_object is None:
            self.abs_tracer_object = AqueousTracer('A520', x, x_aq, y2)
            self.abs_tracer_object.set_normalization_min(abs_prev)
            if self.tracer_view is None:
                self.tracer_view = AqueousTracerView(self, self.abs_tracer_object)
                self.tracer_view.set_data_color(Qt.darkGreen)
        elif self.abs_tracer_object is not None:
            if self.tracer_view is None:
                self.tracer_view = AqueousTracerView(self, self.abs_tracer_object)
                self.tracer_view.set_data_color(Qt.darkGreen)
            elif parent.produced_fluids.cp_or_abs_used:
                self.tracer_view.set_tracer(self.abs_tracer_object)
                self.tracer_view.set_data_color(Qt.darkGreen)

        cg_lyt.addWidget(self.tracer_view, 1, 1, 1, 1)

        button_font = QFont()
        button_font.setPointSize(10)

        if self.viscosity_tracer_object is not None:
            frac_nums_1 = non_nans[0] + 1
            self.v_widget = RetentionWidget(self, self.viscosity_tracer_object, frac_nums_1,
                                            self.viscosity_tracer_object.tracer, 'cP',
                                            parent.produced_fluids.cp_normalized, self.ppm_to_cp_clicked,
                                            lambda: self.calc_from_meas_clicked(3),
                                            self.viscosity_tracer_object.excluded, sf, button_font)
            self.retention_table = self.v_widget.retention_table
            lyt.addWidget(self.v_widget)
            self.v_widget.resident_ppm_edit.setText(str(self.ppm_prev))
            self.tracer_view.staircase_updated.connect(self.retention_table.gray_out_excluded)

        else:
            self.retention_table = None

        if self.abs_tracer_object is not None:
            frac_nums_2 = non_nans2[0] + 1
            self.a_widget = RetentionWidget(self, self.abs_tracer_object, frac_nums_2, self.abs_tracer_object.tracer,
                                            'A520', parent.produced_fluids.abs_normalized, self.ppm_to_abs_clicked,
                                            lambda: self.calc_from_meas_clicked(4), self.abs_tracer_object.excluded,
                                            sf, button_font)
            self.abs_retention_table = self.a_widget.retention_table
            lyt.addWidget(self.a_widget)
            self.a_widget.resident_ppm_edit.setText(str(self.ppm_prev))

            if not v_data_exists or parent.produced_fluids.cp_or_abs_used == 1:
                self.tracer_view.staircase_updated.connect(self.abs_retention_table.gray_out_excluded)

        else:
            self.abs_retention_table = None

        self.ppm_to_cp_tool = None
        self.ppm_to_abs_tool = None
        self.tracer_object = None

        if isinstance(fl, flood.MultiCoreFlood):
            self.va_combo.setCurrentIndex(1)
            self.va_combo.currentIndexChanged.emit()
            self.va_combo.setEnabled(False)
            if fl.effluent.cp_tracer_object is not None:
                self.tracer_view.link_button.setText('Un-link')
        else:
            enable_v = self.viscosity_tracer_object is not None
            enable_a = self.abs_tracer_object is not None
            self.initialize_va_combo(parent.produced_fluids.cp_or_abs_used, enable_v, enable_a)
            if (not parent.produced_fluids.cp_or_abs_used and fl.effluent.cp_tracer_object is not None) or \
                (parent.produced_fluids.cp_or_abs_used and fl.effluent.abs_tracer_object is not None):
                self.tracer_view.link_button.setText('Un-link')

        self.tracer_view.link_tracer.connect(self.link_tracer)

        self.show()

        if not np.isnan(parent.produced_fluids.retention):
            if not parent.produced_fluids.cp_or_abs_used:
                self.v_widget.calc_from_meas_button.clicked.emit()
            else:
                self.a_widget.calc_from_meas_button.clicked.emit()

    def initialize_va_combo(self, cp_or_abs_used: int, enable_v: bool, enable_a: bool):

        self.va_combo.blockSignals(True)
        if not cp_or_abs_used:
            self.v_widget.setVisible(True)
            self.va_combo.setCurrentIndex(0)
            if not enable_a:
                self.va_combo.setEnabled(False)
            else:
                self.a_widget.setVisible(False)
        else:
            self.a_widget.setVisible(True)
            self.va_combo.setCurrentIndex(1)
            if not enable_v:
                self.va_combo.setEnabled(False)
            else:
                self.v_widget.setVisible(False)
        self.va_combo.blockSignals(False)

    def change_palette(self, i):

        if i == 0:
            self.a_widget.setVisible(False)
            self.v_widget.setVisible(True)
            self.tracer_view.staircase_updated.connect(self.retention_table.gray_out_excluded)
            self.tracer_view.staircase_updated.disconnect(self.abs_retention_table.gray_out_excluded)
            self.tracer_view.set_tracer(self.viscosity_tracer_object)
            self.tracer_view.set_data_color(Qt.blue)
            if self.parent().produced_fluids.cp_tracer_object is not None:
                self.tracer_view.link_button.setText('Un-link')
            else:
                self.tracer_view.link_button.setText('Link')

        else:
            self.a_widget.setVisible(True)
            self.v_widget.setVisible(False)
            self.tracer_view.staircase_updated.connect(self.abs_retention_table.gray_out_excluded)
            self.tracer_view.staircase_updated.disconnect(self.retention_table.gray_out_excluded)
            self.tracer_view.set_tracer(self.abs_tracer_object)
            self.tracer_view.set_data_color(Qt.darkGreen)
            if self.parent().produced_fluids.abs_tracer_object is not None:
                self.tracer_view.link_button.setText('Un-link')
            else:
                self.tracer_view.link_button.setText('Link')

    def link_tracer(self, link: bool):

        print('Calling link_tracer({})'.format(link))

        pf = self.parent().produced_fluids
        print(pf, pf.cp_tracer_object)
        source_object = None

        if self.v_widget.isVisible():
            if link:
                source_object = self.viscosity_tracer_object
            pf.cp_tracer_object = source_object

        else:
            if link:
                source_object = self.abs_tracer_object
            pf.abs_tracer_object = source_object
        print(pf, pf.cp_tracer_object, self.viscosity_tracer_object)

    def ppm_to_cp_clicked(self):

        try:
            if self.ppm_to_cp_tool is None:
                self.ppm_to_cp_tool = ToConcentrationTool(self, 'viscosity')
        except Exception as e:
            msg = QMessageBox(parent=self, text=str(e))
            msg.exec_()

    def ppm_to_abs_clicked(self):

        try:
            print('trying ppm to abs', self.ppm_to_abs_tool)
            if self.ppm_to_abs_tool is None:
                self.ppm_to_abs_tool = ToConcentrationTool(self, 'abs')
        except Exception as e:
            msg = QMessageBox(parent=self, text=str(e))
            msg.exec_()

    def calc_from_meas_clicked(self, col: int):

        pf = self.parent().produced_fluids
        if pf.volumes.shape[0] < 3:
            return

        visc = pf.measurements[col, :-1]
        fl = self.parent().parent().flood
        if isinstance(fl, flood.MultiCoreFlood):
            inj_ppm = []
            pln_vol = []
            for rf in fl.ref_floods:
                pln_vol.append(rf.planned_volume)
                ps = rf.fluid.specific_fluid.polymer_solution
                if ps is not None:
                    inj_ppm.append(ps.concentration)
                else:
                    inj_ppm.append(0.)

            inj_ppm_max = max(inj_ppm)
        else:
            inj_ppm = fl.fluid.specific_fluid.polymer_solution.concentration
            inj_ppm_max = inj_ppm
            pln_vol = fl.planned_volume

        if visc.size == 0:
            return
        if col == 3:
            p = copy(pf.ppm_to_cp_fit_params)
            table = self.retention_table
        else:
            p = copy(pf.ppm_to_abs_fit_params)
            table = self.abs_retention_table

        p = np.flipud(p)
        if not p.size:
            return

        if col == 3:
            adj = self.v_widget.meas_adjust_checkbox.checkState()
        else:
            adj = self.a_widget.meas_adjust_checkbox.checkState()

        no_nan = np.where(~np.isnan(visc))
        visc_no_nan = visc[no_nan]
        v_max = np.max(visc_no_nan)

        p_orig = copy(p)
        p[-1] = p_orig[-1] - v_max
        ppm_max = np.real(np.roots(p)[-1])

        if not adj:
            alpha = 1.
        else:
            alpha = ppm_max / inj_ppm_max

        for i in range(np.size(p)):
            p[i] = p[i] * alpha ** (np.size(p) - i - 1)
        p_orig_mod = copy(p)
        p_orig_mod[-1] = p_orig[-1]
        ppm_tracer = np.ones(visc.shape) * np.nan
        for i, v in enumerate(visc):
            if not np.isnan(v):
                p[-1] = p_orig_mod[-1] - v
                new_ppm = np.real(np.roots(p)[-1])
                if new_ppm < 0.:
                    new_ppm = 0.
                ppm_tracer[i] = new_ppm

        table.blockSignals(True)
        for i, val in enumerate(ppm_tracer[no_nan]):
            table.item(i, 2).setText(str(int(val)))
        table.blockSignals(False)

        if col == 3:
            tracer_object = self.viscosity_tracer_object
            color = Qt.blue
        else:
            tracer_object = self.abs_tracer_object
            color = Qt.darkGreen

        old_tracer = tracer_object.tracer.copy()
        tracer_object.tracer = ppm_tracer[no_nan]
        tracer_object.set_normalization_min(self.ppm_prev)
        # self.tracer_view.set_tracer(tracer_object)
        # self.tracer_view.show_data(True)
        # self.tracer_view.set_data_color(color)
        # aq_pv = tracer_object.get_aq_pv()

        # if isinstance(fl, flood.MultiCoreFlood):
        #     inj = 0.
        #     tv = 0.
        #     i = 0
        #     for ippm, plv in zip(inj_ppm, pln_vol):
        #         i += 1
        #         tv += plv
        #         if i < len(inj_ppm) and tv < cv[-1]:
        #             inj += plv * ippm
        #         else:
        #             inj += (cv[-1] - (tv - plv)) * ippm
        # else:
        #     inj = inj_ppm * cv[-1]
        # print('micrograms inj: ', inj)
        pv = fl.get_experiment_total_pore_volume()
        expt = self.parent().parent().parent().experiment
        sw = expt.get_flood_saturation(fl)
        true_aq_pv = pv * sw
        cut = fl.cut_num
        mass = expt.get_core_mass(cut)
        # aq_pv = tracer_object.get_aq_pv(pv * sw)
        rt_pv = tracer_object.get_retention(true_aq_pv)
        if not isinstance(fl, flood.MultiCoreFlood):
            ret = rt_pv * inj_ppm
        else:
            ret = rt_pv * inj_ppm[0]

        # QMessageBox(parent=self, text='Retention is ' + str(np.round(ret / mass)) + ' mcg/g').exec_()

        pf.retention = np.round(ret / mass)
        pf.retention_ppm_min = self.ppm_prev
        ppm = tracer_object.get_interpolated_tracer()
        print('ppm stored in effluent:', ppm)
        pf.ppm = ppm.tolist()
        if col == 3:
            pf.cp_or_abs_used = 0
            pf.cp_normalized = adj
            self.v_widget.set_retention_text(ret / mass, tracer_object.get_d_factor(true_aq_pv))
        else:
            pf.cp_or_abs_used = 1
            pf.abs_normalized = adj
            self.a_widget.set_retention_text(ret / mass, tracer_object.get_d_factor(true_aq_pv))

        tracer_object.tracer = old_tracer
        if col == 3:
            tracer_object.set_normalization_min(self.viscosity_prev)
        else:
            n = p_orig_mod.size
            abs_prev_corrected = np.sum([p_i * self.ppm_prev ** (n - i - 1) for i, p_i in enumerate(p_orig_mod)])
            if self.abs_prev != abs_prev_corrected:
                self.abs_prev = abs_prev_corrected
                print('abs_prev_corrected', abs_prev_corrected)
                tracer_object.set_normalization_min(self.abs_prev)
                self.tracer_view.set_tracer(tracer_object)
                self.tracer_view.set_data_color(color)

    def closeEvent(self, a0: QCloseEvent):
        try:
            self.parent().r_tool = None
            # self.parent().produced_fluids.cp_excluded = self.v_excluded
            if self.abs_retention_table is not None:
                # self.parent().produced_fluids.abs_excluded_row = self.abs_retention_table.excluded_row
                self.parent().produced_fluids.abs_normalized = self.a_widget.meas_adjust_checkbox.isChecked()
            print('cp_is_checked ? {}'.format(self.cp_adjust_checkbox.isChecked()))
            if self.retention_table is not None:
                # self.parent().produced_fluids.cp_excluded_row = self.retention_table.excluded_row
                self.parent().produced_fluids.cp_normalized = self.v_widget.meas_adjust_checkbox.isChecked()
        except Exception as e:
            print(e)


class RetentionWidget(QWidget):

    def __init__(self, parent: RetentionTool, tracer_object: AqueousTracer, frac_nums, meas, units: str,
                 normalized: bool, ppm_to_meas, calc_from_meas, excluded: list, sf: float, button_font: QFont):
        super(RetentionWidget, self).__init__()

        self.tracer_object = tracer_object

        vg_lyt = QGridLayout()
        self.setLayout(vg_lyt)

        self.resident_ppm_widget = QWidget(parent=self)
        self.resident_ppm_widget.setLayout(QHBoxLayout())
        self.resident_ppm_label = QLabel(parent=self.resident_ppm_widget, text='Resident ppm:')
        self.resident_ppm_edit = QLineEdit(parent=self.resident_ppm_widget)
        self.resident_ppm_edit.editingFinished.connect(self.resident_ppm_edit_value_changed)
        self.resident_ppm_widget.layout().addWidget(self.resident_ppm_label)
        self.resident_ppm_widget.layout().addWidget(self.resident_ppm_edit)
        vg_lyt.addWidget(self.resident_ppm_widget, 0, 0, 1, 3)

        self.retention_table = RetentionTable(parent, frac_nums, meas, '[' + units + ']')
        vg_lyt.addWidget(self.retention_table, 1, 0, 1, 3)

        self.meas_adjust_checkbox = QCheckBox(parent=self, text='Normalize?')
        if not normalized:
            self.meas_adjust_checkbox.setChecked(False)
        else:
            self.meas_adjust_checkbox.setChecked(True)
        self.ppm_to_meas_button = QPushButton(parent=self, text='ppm <=> ' + units)
        self.ppm_to_meas_button.setFixedHeight(int(sf * 30))
        self.ppm_to_meas_button.setFont(button_font)
        self.ppm_to_meas_button.clicked.connect(ppm_to_meas)
        self.calc_from_meas_button = QPushButton(parent=self, text='Calc from ' + units)
        self.calc_from_meas_button.setFixedHeight(int(sf * 30))
        self.calc_from_meas_button.setFont(button_font)
        self.calc_from_meas_button.clicked.connect(calc_from_meas)
        vg_lyt.addWidget(self.meas_adjust_checkbox, 2, 0, 1, 1)
        vg_lyt.addWidget(self.ppm_to_meas_button, 2, 1, 1, 1)
        vg_lyt.addWidget(self.calc_from_meas_button, 2, 2, 1, 1)

        self.retention_widget = QWidget(parent=self)
        self.retention_widget.setLayout(QHBoxLayout())
        self.retention_widget.layout().setAlignment(Qt.AlignLeft)
        self.ipv_text = QLabel(parent=self, text='IPV: ')
        self.ipv_text.setFont(button_font)
        self.ipv_text.setFixedWidth(int(sf * 30))
        self.ipv_edit = QLineEdit(parent=self)
        self.ipv_edit.setFont(button_font)
        self.ipv_edit.setFixedWidth(int(sf * 50))
        self.sync_ipv_text_to_stored()
        self.ipv_edit.editingFinished.connect(self.ipv_text_changed)
        self.retention_text = QLabel(parent=self, text='')
        self.retention_text.setFont(button_font)
        self.retention_widget.layout().addWidget(self.ipv_text)
        self.retention_widget.layout().addWidget(self.ipv_edit)
        self.retention_widget.layout().addWidget(self.retention_text)
        vg_lyt.addWidget(self.retention_widget, 3, 0, 1, 3)

        self.meas_adjust_checkbox.setChecked(normalized)
        self.retention_table.gray_out_excluded(excluded)

    def resident_ppm_edit_value_changed(self):

        txt = self.resident_ppm_edit.text()

        try:
            val = float(txt)
            assert val >= 0.
            self.set_ppm_prev(val)

        except ValueError:
            self.sync_ppm_edit()

    def sync_ppm_edit(self):

        self.resident_ppm_edit.setText(str(self.parent().ppm_prev))

    def set_ppm_prev(self, val: float):

        self.parent().ppm_prev = val

    def sync_ipv_text_to_stored(self):

        self.ipv_edit.setText('{:.3f}'.format(self.tracer_object.get_ipv()))

    def ipv_text_changed(self):

        txt = self.ipv_edit.text()

        try:
            val = float(txt)
            self.tracer_object.set_ipv(val)

        except ValueError:
            QMessageBox(parent=self.parent(), text='IPV must be convertible to float.').exec_()

        self.sync_ipv_text_to_stored()

    def set_retention_text(self, ret: float, d: float):

        self.retention_text.setText('  D = {:.3f}, retention = {:.1f} mcg/g'.format(d, ret))


class RetentionTable(QTableWidget):

    def __init__(self, parent: RetentionTool, frac_nums: np.ndarray, data: np.ndarray, header: str):
        super(RetentionTable, self).__init__(parent)

        rc = np.size(data)
        self.frac_nums = frac_nums
        self.data = data
        self.setColumnCount(3)
        self.setRowCount(rc)
        self.setHorizontalHeaderLabels(['#', header, '[ppm]'])
        sf = parent.width() / 900.
        self.setFixedWidth(int(sf * 350))
        font = QFont()
        font.setPointSize(12)
        font.setFamily('Courier')
        self.setFont(font)

        for r in range(rc):
            frac_item = QTableWidgetItem()
            frac_item.setTextAlignment(Qt.AlignCenter)
            visc_item = QTableWidgetItem()
            visc_item.setTextAlignment(Qt.AlignCenter)
            abs_item = QTableWidgetItem()
            abs_item.setTextAlignment(Qt.AlignCenter)
            ppm_item = QTableWidgetItem()
            ppm_item.setTextAlignment(Qt.AlignCenter)
            frac_item.setText(str(frac_nums[r]))
            visc_item.setText(str(data[r]))
            self.setItem(r, 0, frac_item)
            self.setItem(r, 1, visc_item)
            self.setItem(r, 2, ppm_item)

        self.itemClicked.connect(self.is_item_editable)
        self.itemChanged.connect(self.reset_item)

    def reset_item(self, item: QTableWidgetItem):

        if item.column() == 0:
            txt = str(self.frac_nums[item.row()])
        elif item.column() == 1:
            txt = str(self.data[item.row()])
        else:
            txt = ''

        item.setText(txt)

    @staticmethod
    def is_item_editable(item: QTableWidgetItem):
        item.setSelected(False)

    def gray_out_excluded(self, excluded: list):

        if not excluded:
            return

        for row, ex in enumerate(excluded):
            for col in range(self.columnCount()):
                if ex:
                    self.item(row, col).setForeground(Qt.gray)
                else:
                    self.item(row, col).setForeground(Qt.black)


class ToConcentrationTool(QDialog):

    def __init__(self, parent: RetentionTool, mode: str):
        super(ToConcentrationTool, self).__init__(parent=parent)
        self.mode = mode
        sf = parent.width() / 775.
        self.setFixedSize(int(sf * 600), int(sf * 400))
        if mode == 'viscosity':
            prev_order = parent.parent().produced_fluids.ppm_to_cp_fit_params.size - 1
        else:
            prev_order = parent.parent().produced_fluids.ppm_to_abs_fit_params.size - 1

        if prev_order > 0:
            self.fit_order = prev_order
        else:
            self.fit_order = 2

        self.chart = QChart(flags=Qt.WindowFlags())
        self.view = QChartView()
        self.view.setChart(self.chart)
        self.x_axis = QValueAxis()
        self.y_axis = QValueAxis()
        self.chart.setAxisX(self.x_axis)
        self.chart.setAxisY(self.y_axis)

        font = QFont()
        font.setFamily('Courier')
        font.setPointSize(11)
        self.x_axis.setLabelsFont(font)
        self.y_axis.setLabelsFont(font)
        title_font = QFont()
        title_font.setPointSize(12)
        title_font.setBold(True)
        self.x_axis.setTitleText('[ppm]')
        self.x_axis.setTitleFont(title_font)
        self.y_axis.setTitleFont(title_font)
        self.x_axis.setMin(0.)
        self.y_axis.setMin(0.)

        self.data_series = QScatterSeries()
        self.data_series.attachAxis(self.x_axis)
        self.data_series.attachAxis(self.y_axis)
        self.chart.addSeries(self.data_series)

        self.fit_line = QLineSeries()
        self.fit_line.attachAxis(self.x_axis)
        self.fit_line.attachAxis(self.y_axis)
        self.chart.addSeries(self.fit_line)

        self.chart.legend().hide()

        lyt = QGridLayout()
        self.setLayout(lyt)

        lyt.addWidget(self.view, 0, 0, 1, 2)

        self.fit_order_combo = QComboBox(parent=self)
        self.fit_order_combo.addItems(['1', '2', '3'])
        self.fit_order_combo.setCurrentIndex(self.fit_order - 1)
        self.fit_order_combo.currentIndexChanged.connect(self.set_fit_order)
        fit_order_label = QLabel(parent=self, text='Fit Order: ')
        lyt.addWidget(fit_order_label, 1, 0, 1, 1)
        lyt.addWidget(self.fit_order_combo, 1, 1, 1, 1)

        if mode == 'viscosity':
            self.setWindowTitle('[cP] to [ppm] Tool')
            self.y_axis.setTitleText('[cP]')
            ppm_to_m_data = self.parent().parent().produced_fluids.ppm_to_cp
        else:
            self.setWindowTitle('A520 to [ppm] Tool')
            self.y_axis.setTitleText('A520')
            ppm_to_m_data = self.parent().parent().produced_fluids.ppm_to_abs
        if ppm_to_m_data.size > 0:
            ppm_data = ppm_to_m_data[0, :]
            m_data = ppm_to_m_data[1, :]
        else:
            ppm_data = np.array([])
            m_data = np.array([])

        self.to_ppm_table = ToConcentrationTable(self, ppm_data[:], m_data[:])
        if ppm_to_m_data.size > 0:
            self.update_plot()

        lyt.addWidget(self.to_ppm_table, 0, 2, 1, 1)

        self.import_from_rheology_button = QPushButton(parent=self, text='Import from Rheology')
        self.import_from_rheology_button.clicked.connect(self.import_from_rheology)
        lyt.addWidget(self.import_from_rheology_button, 1, 2, 1, 1)

        self.show()

    def set_fit_order(self, i: int):

        self.fit_order = i + 1
        self.update_plot()

    def import_from_rheology(self):

        try:
            fl = self.parent().parent().parent().flood
            fluids = fl.get_fluids_list()
            polymers = []
            brines = []
            for fld in fluids:
                if fld.get_fluid_type() == FluidType.POLYMER_SOLUTION:
                    polymers.append(fld.specific_fluid.polymer_solution.primary_additive)
                    brines.append(fld.specific_fluid.polymer_solution.base_fluid)

            if not polymers:
                return

            p = polymers[0]

            if not all(utils.list_element_wise_compare(polymers, p)):
                return

            pf = self.parent().parent().produced_fluids
            t = pf.temperature
            shear = pf.gamma

            cs = []
            vs = []

            for b in brines:
                rs = p.get_matching_rheologies(b)
                if rs:
                    cs, vs = p.get_viscosity_vs_concentration_data_interp(b, t, shear)
                    if cs:
                        cs = [0., *cs]
                        vs = [b.get_viscosity(t), *vs]
                        break

            if not cs:
                return

            self.to_ppm_table.ppm_data = np.array(cs)
            self.to_ppm_table.m_data = np.array(vs)
            self.to_ppm_table.setRowCount(1)
            self.to_ppm_table.set_from_data()
            self.update_plot()

        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()

    def update_plot(self):

        if self.data_series.pointsVector():
            self.data_series.removePoints(0, len(self.data_series.pointsVector()))
        if self.fit_line.pointsVector():
            self.fit_line.removePoints(0, len(self.fit_line.pointsVector()))
        self.data_series.attachAxis(self.x_axis)
        self.data_series.attachAxis(self.y_axis)
        self.fit_line.attachAxis(self.x_axis)
        self.fit_line.attachAxis(self.y_axis)

        ppm = self.to_ppm_table.ppm_data
        meas = self.to_ppm_table.m_data
        non_nan = np.logical_and(~np.isnan(ppm), ~np.isnan(meas))
        ppm = ppm[non_nan]
        meas = meas[non_nan]
        if self.mode == 'viscosity':
            self.parent().parent().produced_fluids.ppm_to_cp = np.array([ppm.tolist(), meas.tolist()])
        else:
            self.parent().parent().produced_fluids.ppm_to_abs = np.array([ppm.tolist(), meas.tolist()])

        self.data_series.append(utils.series_to_polyline(ppm.tolist(), meas.tolist()))
        self.data_series.setMarkerSize(12.)
        self.data_series.setColor(Qt.black)
        if self.mode == 'viscosity':
            self.parent().parent().produced_fluids.ppm_to_cp_fit_params = np.array([])
        else:
            self.parent().parent().produced_fluids.ppm_to_abs_fit_params = np.array([])

        if ppm.tolist():
            self.x_axis.setMax(np.max(ppm))
            self.y_axis.setMax(np.max(meas))
            used_order = min(self.fit_order, len(ppm.tolist()) - 1)
            if used_order == 1:
                z, _ = curve_fit(ToConcentrationTool.general_polynomial, xdata=ppm, ydata=meas,
                              p0=[1., 0.01], bounds=(0., np.inf))
            elif used_order == 2:
                z, _ = curve_fit(ToConcentrationTool.general_polynomial, xdata=ppm, ydata=meas,
                              p0=[1., 0.01, 1.e-6], bounds=(0., np.inf))
            elif used_order == 3:
                z, _ = curve_fit(ToConcentrationTool.general_polynomial, xdata=ppm, ydata=meas,
                              p0=[1., 0.01, 1.e-6, 1.e-9], bounds=(0., np.inf))
            else:
                return

            if self.mode == 'viscosity':
                self.parent().parent().produced_fluids.ppm_to_cp_fit_params = z
            else:
                self.parent().parent().produced_fluids.ppm_to_abs_fit_params = z
            fit_ppm = np.linspace(start=ppm[0], stop=ppm[-1])
            fit_meas = ToConcentrationTool.general_polynomial(fit_ppm, *z)
            self.fit_line.append(utils.series_to_polyline(fit_ppm.tolist(), fit_meas.tolist()))
            self.fit_line.setColor(Qt.red)

    @staticmethod
    def general_polynomial(x, *args):
        ans = 0.
        for i, arg in enumerate(args):
            ans += arg * np.power(x, i)
        return ans

    def closeEvent(self, a0: QCloseEvent):
        if self.mode == 'viscosity':
            self.parent().ppm_to_cp_tool = None
        else:
            self.parent().ppm_to_abs_tool = None


class ToConcentrationTable(QTableWidget):

    rc_changed = pyqtSignal(int)
    data_changed = pyqtSignal()

    def __init__(self, parent: ToConcentrationTool, ppm_data: np.array, m_data: np.array):
        super(ToConcentrationTable, self).__init__(parent=parent)
        self.setRowCount(1)
        self.setColumnCount(2)
        sf = parent.width() / 600.
        self.setFixedWidth(int(sf * 200))
        for i in range(2):
            self.setColumnWidth(i, int(sf * 90))
        if parent.mode == 'viscosity':
            self.setHorizontalHeaderLabels(['[ppm]', '[cP]'])
        else:
            self.setHorizontalHeaderLabels(['[ppm]', 'A520'])
        self.ppm_data = ppm_data
        self.m_data = m_data

        self.itemChanged.connect(self.check_item)
        self.rc_changed.connect(self.update_rows)
        self.data_changed.connect(parent.update_plot)

        if ppm_data.size != 0:
            self.set_from_data()

    def set_from_data(self):

        for i in range(self.ppm_data.size):
            for j in range(2):
                print(i, j)
                item = QTableWidgetItem()
                if j == 0:
                    item.setText(str(self.ppm_data[i]))
                else:
                    item.setText(str(self.m_data[i]))
                try:
                    self.blockSignals(True)
                    self.setItem(i, j, item)
                    self.blockSignals(False)
                except Exception as e:
                    print(e)
            self.blockSignals(True)
            self.setRowCount(self.rowCount() + 1)
            self.blockSignals(False)
        self.ppm_data = np.append(self.ppm_data, np.nan)
        self.m_data = np.append(self.m_data, np.nan)

    def check_item(self, item: QTableWidgetItem):

        try:
            val = float(item.text())
            assert val >= 0.
        except Exception as e:
            print(e)
            if item.row() == self.rowCount() - 1:
                del item
            else:
                if item.column() == 0:
                    if np.isnan(self.ppm_data[item.row()]):
                        item.setText('')
                    else:
                        item.setText(str(self.ppm_data[item.row()]))
                else:
                    if np.isnan(self.m_data[item.row()]):
                        item.setText('')
                    else:
                        item.setText(str(self.m_data[item.row()]))
            return

        if item.row() == self.rowCount() - 1:
            self.rc_changed.emit(self.rowCount() + 1)

        if item.column() == 0:
            self.ppm_data[item.row()] = val
        else:
            self.m_data[item.row()] = val

        self.data_changed.emit()

    def update_rows(self, nrc: int):

        rc = self.rowCount()

        if nrc < rc:
            self.ppm_data = self.ppm_data[:nrc]
            self.m_data = self.m_data[:nrc]
        else:
            for i in range(rc, nrc):
                item0 = self.item(i - 1, 0)
                item1 = self.item(i - 1, 1)

                if item0 is not None:
                    ppm = np.append(self.ppm_data, float(item0.text()))
                else:
                    ppm = np.append(self.ppm_data, np.nan)
                self.ppm_data = ppm

                if item1 is not None:
                    meas = np.append(self.m_data, float(self.item1.text()))
                else:
                    meas = np.append(self.m_data, np.nan)
                self.m_data = meas

        self.setRowCount(nrc)

    def keyPressEvent(self, e: QKeyEvent):

        if e.key() == Qt.Key_Delete:
            item = self.currentItem()
            if item is None:
                return

            r = item.row()

            if r == self.rowCount() - 1:
                return

            if item.column() == 0:
                self.ppm_data[r] = np.nan
            else:
                self.m_data[r] = np.nan

            item = self.takeItem(r, item.column())
            del item

            if np.isnan(self.ppm_data[r]) and np.isnan(self.m_data[r]):
                nr = self.rowCount()
                self.removeRow(r)
                if r != nr - 1:
                    pd = np.append(self.ppm_data[:r], self.ppm_data[r + 1:])
                    md = np.append(self.m_data[:r], self.m_data[r + 1:])
                else:
                    pd = self.ppm_data[:r]
                    md = self.m_data[:r]
                self.ppm_data = pd
                self.m_data = md

            self.data_changed.emit()

        else:
            super(ToConcentrationTable, self).keyPressEvent(e)


class AlkaliConsumptionTool(QDialog):

    def __init__(self, parent: ProducedFluidsView):
        super(AlkaliConsumptionTool, self).__init__(parent=parent)

        self.setWindowTitle('Alkali Consumption Tool')
        layout = QVBoxLayout()
        h_layout = QHBoxLayout()
        h_layout_2 = QHBoxLayout()
        self.setLayout(layout)
        self.consumption_text = QLabel(parent=self, text='Moles OH- Recovered: ')
        font = self.consumption_text.font()
        font.setPointSize(12)
        self.consumption_text.setFont(font)
        h_layout.addWidget(self.consumption_text)
        self.send_consumption_button = QPushButton(parent=self, text='Set as alkali consumption')
        self.send_consumption_button.setFont(font)
        self.send_consumption_button.clicked.connect(self.update_flood_alkali_consumption)

        h_layout_2.addWidget(self.send_consumption_button)

        layout.addLayout(h_layout)
        layout.addLayout(h_layout_2)

        sf = parent.width() / 1010.
        self.setFixedSize(int(sf * 550), int(sf * 550))

        self.cv = parent.produced_fluids.cumulative_volume()
        self.cv_aq = parent.produced_fluids.cumulative_aq_volume()
        name = '[OH-]'

        if parent.parent().flood.alkali_consumption_tracer_object is None:

            ph = parent.produced_fluids.measurements[0, :]
            tracer = np.power(10., ph - 14.)

            non_nan = np.where(~np.isnan(tracer))

            if not non_nan[0].tolist():
                QMessageBox(parent=parent, text='No non-nan data available.').exec_()
                return

            self.tracer_object = AqueousTracer(name, self.cv.copy(), self.cv_aq.copy(), tracer.copy())
            self.tracer_object.set_normalization_max(1.)
            self.send_consumption_button.setEnabled(False)

        else:

            self.tracer_object = parent.parent().flood.alkali_consumption_tracer_object

        self.tracer_view = AqueousTracerView(parent=self, tracer=self.tracer_object)
        self.tracer_view.y_axis.setMax(np.max(self.tracer_object.tracer))
        self.tracer_view.y_axis.setTitleText(name + ' [M]')

        self.tracer_view.staircase_updated.connect(self.update_consumption)
        self.tracer_view.link_tracer.connect(self.link_button_clicked)

        layout.addWidget(self.tracer_view)

        self.consumption = np.nan
        self.update_consumption([])
        self.show()

    def update_consumption_text(self):

        inj = self.calc_injected_recovered()
        txt = 'Moles OH- Inj.: {:.2E}, Rec.: {:.2E}, Consumption: {:.1f}%'
        self.consumption_text.setText(txt.format(inj, self.consumption, 100. * (inj - self.consumption) / inj))

    def calc_consumption(self):

        return 0.001 * self.tracer_object.get_area_above_below_tracer(False)

    def calc_injected_recovered(self):

        fl = self.parent().parent().flood
        sw = fl.experiment.get_flood_saturation(fl)
        pv = fl.get_experiment_total_pore_volume()
        f_list = fl.get_fluids_list()
        volume_injected = self.tracer_object.get_inc_cv()[-1]
        volume_remaining = pv * sw
        injected_volume_recovered = volume_injected - volume_remaining

        inj = 0.
        tot_vol = 0.

        if isinstance(fl, flood.MultiCoreFlood):
            fl_list = fl.ref_floods
            for fli, fi in zip(fl_list[:-1], f_list[:-1]):
                if tot_vol + fli.planned_volume <= injected_volume_recovered:
                    inj += 0.001 * fli.planned_volume * np.power(10., fi.specific_fluid.pH - 14.)
                    tot_vol += fli.planned_volume
                else:
                    inj += 0.001 * (injected_volume_recovered - tot_vol) * np.power(10., fi.specific_fluid.pH - 14.)
                    tot_vol = injected_volume_recovered
                    break

        inj += 0.001 * (injected_volume_recovered - tot_vol) * np.power(10., f_list[-1].specific_fluid.pH - 14.)

        return inj

    def link_button_clicked(self, link: bool):

        if link:
            self.set_flood_tracer_object()
            self.send_consumption_button.setEnabled(True)
        else:
            self.clear_flood_tracer_object()
            self.send_consumption_button.setEnabled(False)

    def update_flood_alkali_consumption(self):

        self.parent().parent().flood.alkali_injected_recovered = self.calc_injected_recovered()
        self.parent().parent().flood.alkali_consumption = self.consumption

    def set_flood_tracer_object(self):

        self.parent().parent().flood.alkali_consumption_tracer_object = self.tracer_object

    def clear_flood_tracer_object(self):

        self.parent().parent().flood.alkali_consumption_tracer_object = None

    def update_consumption(self, _):

        self.consumption = self.calc_consumption()
        self.update_consumption_text()

    def closeEvent(self, _):

        try:
            self.parent().c_tool = None
        except Exception as e:
            print(e)


if __name__ == '__main__':

    from PyQt5.Qt import QApplication
    import sys

    app = QApplication(sys.argv)

    produced_fluids = ProducedFluids(True, True)
    # produced_fluids.volumes = np.array([[5., 4.], [10., 10.], [5., 6.]])
    produced_fluids_view = ProducedFluidsView(produced_fluids)

    sys.exit(app.exec_())
