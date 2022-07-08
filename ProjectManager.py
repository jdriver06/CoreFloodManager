
from PyQt5.QtWidgets import QDialog, QWidget, QListWidget, QHBoxLayout, QVBoxLayout, QLabel, QPushButton, \
    QGroupBox, QLineEdit, QGridLayout, QComboBox, QMessageBox, QListWidgetItem, QMainWindow
from PyQt5.QtGui import QCloseEvent, QMouseEvent
from PyQt5.Qt import QIcon, pyqtSignal, QApplication
from PyQt5.QtCore import QModelIndex
import j_utils as utils
from core import Core, CoreSection, CoreView, CoreViewForm, CoreEntryForm, CoreHolder, CoreSaveDlg, \
    CoreSaveStructure, JQTextEdit
import fluid
import flood
import numpy as np
import flood_experiment
import pickle as pkl
from functools import partial
import os.path as osp
from os import remove
from shutil import copyfile
from filelock import FileLock
from copy import deepcopy


__version__ = '0.10.0'

# Patched issue with extra nan entry in oil if temps by discarding in polymer solution brine_if.

version_list = ['0.1.0',
                '0.1.1',
                '0.2.0',
                '0.2.1',
                '0.3.0',
                '0.3.1',
                '0.4.0',
                '0.5.0',
                '0.5.1',
                '0.6.0',
                '0.6.1',
                '0.6.2',
                '0.6.3',
                '0.7.0',
                '0.7.1',
                '0.8.0',
                '0.8.1',
                '0.9.0',
                '0.9.1',
                '0.9.2',
                '0.9.3',
                '0.9.4',
                '0.9.5',
                '0.9.6',
                '0.9.7',
                '0.10.0']


class ProjectManager:

    def __init__(self):
        self.core_pickles = []
        self.name = ''
        self.version = __version__


class NameDialog(QDialog):

    def __init__(self, parent, fields: list, data_types: list):
        super(NameDialog, self).__init__(parent=parent)

        if not fields:
            self.close()
            return

        if len(fields) != len(data_types) + 1:
            self.close()
            return

        icon = QIcon('Coreholder Cropped.jpg')
        self.setWindowIcon(icon)
        self.setWindowTitle(' ')

        lyt = QVBoxLayout()
        self.setLayout(lyt)

        self.labels = []
        self.edits = []
        self.data_types = data_types

        for field in fields:
            if isinstance(field, list):
                self.labels.append(QLabel(parent=self, text=field[0] + ':'))
                self.edits.append(QComboBox(parent=self))
                self.edits[-1].addItems(field[1:])
            else:
                self.labels.append(QLabel(parent=self, text=field + ':'))
                self.edits.append(QLineEdit(parent=self))

            lyt.addWidget(self.labels[-1])
            lyt.addWidget(self.edits[-1])

        self.b_layout = QHBoxLayout(self)
        lyt.addLayout(self.b_layout)

        submit = QPushButton(parent=self, text='Submit')
        submit.clicked.connect(self.submit)

        self.b_layout.addWidget(submit)

        self.show()

    def submit(self):
        name = self.edits[0].text()
        if not name:
            return
        else:
            args = []
            if self.data_types:
                for i, data_type in enumerate(self.data_types):
                    e = self.edits[i+1]
                    try:
                        if isinstance(e, QLineEdit):
                            args.append(data_type(e.text()))
                        elif isinstance(e, QComboBox):
                            args.append(e.itemText(e.currentIndex()))
                    except Exception as e:
                        print(e)
                        return

            self.parent().submit_name(name, *args)
            self.close()

    def closeEvent(self, a0: QCloseEvent):
        p = self.parent()
        if not p.isVisible():
            p.close()


class CFM_NameDialog(NameDialog):

    def __init__(self, parent):
        super(CFM_NameDialog, self).__init__(parent, ['Project Name'], [])

        load = QPushButton(parent=self, text='Load')
        load.clicked.connect(self.load)
        self.b_layout.addWidget(load)

    def load(self):

        file_name = utils.file_open_dlg('', '(Core Flood Manager) *.cfm', False)
        if not file_name or osp.isfile(file_name + '.lock'):
            return
        else:
            try:
                lock = FileLock(file_name + '.lock')
                lock.acquire()
                pickle_in = open(file_name, 'rb')
            except Exception as e:
                msg = QMessageBox(parent=self, title='File load error',
                                  text=e)
                msg.exec_()
                return
            try:
                obj = pkl.load(pickle_in)
                self.parent().load_from_pickle(obj, pickle_in, lock)
            except Exception as e:
                print(e)
            self.parent().setVisible(True)
            self.close()


class ProjectManagerView(QDialog):

    def __init__(self):
        super(ProjectManagerView, self).__init__()
        self.source_file = None
        self.file_lock = None

        CFM_NameDialog(self)

        self.project_manager = ProjectManager()

        icon = QIcon('Coreholder Cropped.jpg')
        self.setWindowIcon(icon)
        sf = 92. / QApplication.primaryScreen().physicalDotsPerInchY()
        orientation = QApplication.primaryScreen().primaryOrientation()
        avail_size = QApplication.primaryScreen().availableSize()
        physical_size = QApplication.primaryScreen().physicalSize()
        print(orientation, avail_size, physical_size)
        self.sf = sf
        self.setFixedSize(int(sf * 900), int(sf * 700))

        layout = QGridLayout()
        self.setLayout(layout)

        self.plug_list = utils.SignalListManager()
        self.plug_list.add_signal_list()
        self.plug_list_widget = utils.SignalListManagerWidget(self, self.plug_list, ['Plugs'], [[NameDialog], [CoreSection]],
                                                         ['Name', 'Diameter [cm]', 'Length [cm]', 'Initial Mass [g]',
                                                          ['Lithology', 'sandstone', 'limestone', 'dolomite'],
                                                          chr(981) + ' [frac] (est.)'],
                                                         [float, float, float, str, float])
        self.oil_list = utils.SignalListManager()
        self.oil_list.add_signal_list(3)
        self.oil_list_widget = utils.SignalListManagerWidget(self, self.oil_list,
                                                             ['Oils', 'Samples/Surrogates', 'Diluents'],
                                                        [[NameDialog, NameDialog, NameDialog],
                                                         [fluid.Oil, fluid.OilSample, fluid.Diluent]], ['Name'], [])
        self.oil_list_widget.pm_list.signal_lists[1].set_ref_lists([self.oil_list_widget.pm_list.signal_lists[0]])

        self.brine_list = utils.SignalListManager()
        self.brine_list.add_signal_list()
        self.brine_list_widget = utils.SignalListManagerWidget(self, self.brine_list, ['Brines'],
                                                               [[fluid.BrineTool], [fluid.Brine]], ['Name'], [])

        self.chemical_list = utils.SignalListManager()
        self.chemical_list.add_signal_list(6)
        self.chemical_list_widget = utils.SignalListManagerWidget(self, self.chemical_list, ['Polymers', 'Co-solvents',
                                                                                             'Surfactants',
                                                                                             'Scavengers', 'Other',
                                                                                             'Formulations'],
                                                             [[NameDialog, NameDialog, NameDialog, NameDialog,
                                                               NameDialog, fluid.BrineInjectionFluidFormulationViewWrapper],
                                                              [fluid.Polymer, fluid.CoSolvent, fluid.Surfactant,
                                                               fluid.Scavenger, fluid.Other, fluid.Formulation]],
                                                             ['Name'], [])
        self.chemical_list.signal_lists[5].set_ref_lists([self.brine_list.signal_lists[0]])

        self.core_list = utils.SignalListManager()
        self.core_list.add_signal_list()
        self.core_list_widget = utils.SignalListManagerWidget(self, self.core_list, ['Cores'], [[NameDialog], [Core]],
                                                         ['Name'], [])
        self.core_list_widget.pm_list.signal_lists[0].set_ref_lists([self.plug_list_widget.pm_list.signal_lists[0]])

        self.flood_experiment_list = utils.SignalListManager()
        self.flood_experiment_list.add_signal_list()
        self.flood_experiment_list_widget = utils.SignalListManagerWidget(self, self.flood_experiment_list, ['Experiments'],
                                                                     [[NameDialog],
                                                                      [flood_experiment.FloodExperiment]], ['Name'], [])
        self.flood_experiment_list_widget.pm_list.signal_lists[0].\
            set_ref_lists([self.core_list_widget.pm_list.signal_lists[0]])

        self.gas_list = utils.SignalListManager()
        self.gas_list.add_signal_list()
        self.gas_list_widget = utils.SignalListManagerWidget(self, self.gas_list, ['Gases'], [[NameDialog], [fluid.Gas]],
                                                        ['Name'], [])

        layout.addWidget(self.plug_list_widget, 0, 0, 1, 1)
        layout.addWidget(self.core_list_widget, 0, 1, 1, 1)
        layout.addWidget(self.flood_experiment_list_widget, 0, 2, 1, 1)

        layout.addWidget(self.gas_list_widget, 1, 0, 1, 1)
        layout.addWidget(self.oil_list_widget, 1, 1, 1, 2)

        layout.addWidget(self.brine_list_widget, 2, 0, 1, 1)
        layout.addWidget(self.chemical_list_widget, 2, 1, 1, 2)

        self.save_button = QPushButton(parent=self, text='Save')
        self.save_button.clicked.connect(lambda: self.dump_pickle(True))
        layout.addWidget(self.save_button, 3, 0, 1, 3)

    def submit_name(self, name: str):
        self.project_manager.name = name
        self.setWindowTitle('Project Manager v' + __version__ + ': ' + self.project_manager.name)
        self.show()

    def create_pickle(self):
        return ProjectManagerSaveStructure(self)

    def load_from_pickle(self, pickle, source_file, file_lock):

        self.source_file = source_file
        self.file_lock = file_lock
        self.project_manager = pickle.project_manager
        self.submit_name(self.project_manager.name)
        self.plug_list_widget.clear_items()
        self.plug_list = pickle.plug_list
        self.plug_list_widget.pm_list = pickle.plug_list
        self.plug_list_widget.populate_from_pm_list()
        for c_list in self.plug_list_widget.pm_list.signal_lists:
            for obj in c_list.objects:
                obj.project_manager_list = self.plug_list_widget

        self.oil_list_widget.clear_items()
        self.oil_list = pickle.oil_list
        # self.oil_list.signal_lists[1].set_ref_lists([self.oil_list.signal_lists[0],
        #                                              self.oil_list.signal_lists[2]])
        self.oil_list_widget.pm_list = pickle.oil_list
        self.oil_list_widget.populate_from_pm_list()
        for c_list in self.oil_list_widget.pm_list.signal_lists:
            for obj in c_list.objects:
                obj.project_manager_list = self.oil_list_widget

        self.flood_experiment_list_widget.clear_items()
        self.flood_experiment_list = pickle.flood_experiment_list
        self.flood_experiment_list_widget.pm_list = pickle.flood_experiment_list
        self.flood_experiment_list_widget.populate_from_pm_list()
        # import numpy as np
        for c_list in self.flood_experiment_list_widget.pm_list.signal_lists:
            for obj in c_list.objects:
                obj.project_manager_list = self.flood_experiment_list_widget

        self.brine_list_widget.clear_items()
        self.brine_list = pickle.brine_list
        self.brine_list_widget.pm_list = pickle.brine_list
        self.brine_list_widget.populate_from_pm_list()
        for c_list in self.brine_list_widget.pm_list.signal_lists:
            for obj in c_list.objects:
                obj.project_manager_list = self.brine_list_widget

        self.core_list_widget.clear_items()
        self.core_list = pickle.core_list
        self.core_list_widget.pm_list = pickle.core_list
        self.core_list_widget.populate_from_pm_list()
        for c_list in self.core_list_widget.pm_list.signal_lists:
            for obj in c_list.objects:
                obj.project_manager_list = self.core_list_widget

        self.chemical_list_widget.clear_items()
        self.chemical_list = pickle.chemical_list
        if len(self.chemical_list.signal_lists) == 5:
            self.chemical_list.add_signal_list()
        self.chemical_list_widget.pm_list = self.chemical_list
        self.chemical_list.signal_lists[5].set_ref_lists([self.brine_list.signal_lists[0]])
        self.chemical_list_widget.populate_from_pm_list()
        for c_list in self.chemical_list.signal_lists:
            for obj in c_list.objects:
                if obj is not None:
                    obj.project_manager_list = self.chemical_list_widget

        self.gas_list_widget.clear_items()
        self.gas_list = pickle.gas_list
        self.gas_list_widget.pm_list = pickle.gas_list
        self.gas_list_widget.populate_from_pm_list()
        for c_list in self.gas_list_widget.pm_list.signal_lists:
            for obj in c_list.objects:
                obj.project_manager_list = self.gas_list_widget

        if self.project_manager.version > '0.4.0':
            self.project_manager.version = '0.4.0'
        self.upgrade_version_to_current()
        # self.project_manager.version = '0.4.0'

    def upgrade_version_to_current(self):

        print(self.project_manager.version)

        if not hasattr(self.project_manager, 'version'):
            self.project_manager.version = '0.1.0'
        ind = version_list.index(str(self.project_manager.version))
        print(version_list)
        for i in range(ind, len(version_list) - 1):
            self.version_increment(i)

    def version_increment(self, i: int):
        self.project_manager.version = version_list[i + 1]
        print(version_list[i])

        if version_list[i] == '0.1.0':
            pass

        if version_list[i] == '0.1.1':
            for obj in self.flood_experiment_list.signal_lists[0].objects:
                for fl in obj.floods:
                    fl.effluent.ppm_to_cp = np.array([[], []])
                    fl.effluent.ppm_to_abs = np.array([[], []])
                    fl.effluent.cp_excluded_row = -1
                    fl.effluent.abs_excluded_row = -1
                    fl.effluent.cp_excluded = []
                    fl.effluent.abs_excluded = []
                    fl.effluent.cp_normalized = 0
                    fl.effluent.abs_normalized = 0
                    fl.effluent.cp_or_abs_used = 0
                    fl.effluent.ppm = np.array([])
                    fl.effluent.retention = np.nan
                    fl.effluent.ppm_to_cp_fit_params = np.array([])
                    fl.effluent.ppm_to_abs_fit_params = np.array([])
                    m = deepcopy(fl.effluent.measurements)
                    a = np.ones((1, m.shape[1]))
                    a[:] = np.nan
                    a = a[0].tolist()
                    m = m.tolist()
                    m.append(a)
                    fl.effluent.measurements = np.array(m)

        if version_list[i] == '0.2.0':
            pass

        if version_list[i] == '0.2.1':
            for obj in self.chemical_list.signal_lists[0].objects:
                if not hasattr(obj, 'rheologies_list'):
                    obj.rheologies_list = utils.SignalListManager()
                    obj.rheologies_list.add_signal_list()
                if not hasattr(obj, 'base_fluids_list'):
                    obj.base_fluids_list = utils.SignalListManager()
                    obj.base_fluids_list.add_signal_list()
            for obj in self.flood_experiment_list.signal_lists[0].objects:
                for fl in obj.floods:
                    if not hasattr(fl, 'krw'):
                        fl.krw = np.nan
                    if isinstance(fl.fluid, fluid.Brine):
                        fl.fluid = fluid.InjectionFluid('', [fluid.Brine('')], [[], 'brine-based'])
                    for attr in fl.effluent.__dict__.keys():
                        print(attr)
            for i, obj in enumerate(self.brine_list.signal_lists[0].objects):
                if not hasattr(obj, 'add_salt_composition'):
                    new_brine = fluid.Brine(obj.name)
                    new_brine.project_manager_list = obj.project_manager_list
                    new_brine.ref_viscosity = obj.ref_viscosity
                    for ion, value in obj.composition.items():
                        new_brine.composition[ion] = value
                    for attr in new_brine.__dict__.keys():
                        obj.__dict__[attr] = new_brine.__dict__[attr]
                else:
                    print(obj.name, obj.add_salt_composition)

        if version_list[i] == '0.3.0':
            for obj in self.flood_experiment_list.signal_lists[0].objects:
                for fl in obj.floods:
                    if not hasattr(fl, 'y_lims'):
                        fl.y_lims = [np.nan, np.nan]
                    vols = fl.effluent.volumes
                    fl.effluent.volumes[0, :] = np.round(vols[0, :], 3)
                    fl.effluent.volumes[1, :] = np.round(vols[1, :], 3)
                    fl.effluent.volumes[2, :] = np.round(vols[1, :] - vols[0, :], 3)

        if version_list[i] == '0.3.1':
            for obj in self.flood_experiment_list.signal_lists[0].objects:
                if not hasattr(obj, 'floods_list'):
                    obj.floods_list = utils.SignalListManager()
                    obj.floods_list.add_signal_list()
                if not hasattr(obj, 'multi_floods_list'):
                    obj.multi_floods_list = utils.SignalListManager()
                    obj.multi_floods_list.add_signal_list()
                    obj.multi_floods_list.signal_lists[0].set_ref_lists([obj.floods_list.signal_lists[0]])
                for mf in obj.multi_floods_list.signal_lists[0].objects:
                    mf.n_phases = mf.ref_floods[0].n_phases
                    if not hasattr(mf, 'excluded'):
                        mf.excluded = [[], [], [], [], []]
                    if not mf.excluded:
                        mf.excluded = [[], [], [], [], []]
                for fl in obj.floods:
                    if not hasattr(fl, 'excluded'):
                        fl.excluded = [[], [], [], [], []]
                    elif not fl.excluded:
                        fl.excluded = [[], [], [], [], []]
                    if not hasattr(fl, 'experiment'):
                        fl.experiment = obj

        if version_list[i] == '0.4.0':
            for obj in self.flood_experiment_list.signal_lists[0].objects:
                for fld in obj.floods:
                    if not hasattr(fld, 'rf_source_flood'):
                        fld.rf_source_flood = fld
                    if not hasattr(fld, 'rf_referencer_floods'):
                        fld.rf_referencer_floods = []
                    if not hasattr(fld, 'permeability_viscosity'):
                        fld.permeability_viscosity = np.nan
                    if not hasattr(fld, 'y_lims_rf'):
                        fld.y_lims_rf = [np.nan, np.nan]
                    if not hasattr(fld, 'plateau_whole_rf'):
                        fld.plateau_whole_rf = []
                    if not hasattr(fld, 'plateau_sec1_rf'):
                        fld.plateau_sec1_rf = []
                    if not hasattr(fld, 'plateau_sec2_rf'):
                        fld.plateau_sec2_rf = []
                    if not hasattr(fld, 'plateau_sec3_rf'):
                        fld.plateau_sec3_rf = []
                    if not hasattr(fld, 'plateau_sec4_rf'):
                        fld.plateau_sec4_rf = []
                    if not hasattr(fld, 'rf_data'):
                        fld.rf_data = [np.ndarray([]), np.ndarray([]), np.ndarray([]), np.ndarray([]), np.ndarray([])]
                for mf in obj.multi_floods_list.signal_lists[0].objects:
                    if not hasattr(mf, 'rf_source_flood'):
                        mf.rf_source_flood = mf
                    if not hasattr(mf, 'rf_referencer_floods'):
                        mf.rf_referencer_floods = []
                    if not hasattr(mf, 'permeability_viscosity'):
                        mf.permeability_viscosity = np.nan
                    if not hasattr(mf, 'y_lims_rf'):
                        mf.y_lims_rf = [np.nan, np.nan]
                    if not hasattr(mf, 'plateau_whole_rf'):
                        mf.plateau_whole_rf = []
                    if not hasattr(mf, 'plateau_sec1_rf'):
                        mf.plateau_sec1_rf = []
                    if not hasattr(mf, 'plateau_sec2_rf'):
                        mf.plateau_sec2_rf = []
                    if not hasattr(mf, 'plateau_sec3_rf'):
                        mf.plateau_sec3_rf = []
                    if not hasattr(mf, 'plateau_sec4_rf'):
                        mf.plateau_sec4_rf = []
                    if not hasattr(mf, 'rf_data'):
                        mf.rf_data = [np.ndarray([]), np.ndarray([]), np.ndarray([]), np.ndarray([]), np.ndarray([])]

        if version_list[i] == '0.5.0':
            pass

        if version_list[i] == '0.5.1':
            for expt in self.flood_experiment_list.signal_lists[0].objects:
                for fld in expt.floods:
                    vols = fld.effluent.volumes
                    s = vols.shape
                    if s[0] == 3:
                        new_vols = np.zeros((4, s[1]))
                        new_vols[1, :] = np.nan
                        new_vols[0, :] = vols[0, :]
                        new_vols[2:, :] = vols[1:, :]
                        fld.effluent.volumes = new_vols
                    if not hasattr(fld.effluent, 'swelling_factor'):
                        fld.effluent.swelling_factor = 1.
                    if not hasattr(fld.effluent, 'revisions'):
                        fld.effluent.revisions = np.zeros((4, s[1]))
                        fld.effluent.revisions[1, :] = np.nan
                    if not hasattr(fld.effluent, 'micro_emulsions'):
                        fld.effluent.micro_emulsions = np.zeros((2, s[1]))
                        fld.effluent.micro_emulsions[:, :] = np.nan
                    if not hasattr(fld.effluent, 'ions') or not fld.effluent.ions.size:
                        fld.effluent.ions = np.zeros((10, s[1]))
                        fld.effluent.ions[:, :] = np.nan
                    if not hasattr(fld.effluent, 'dead_volume'):
                        fld.effluent.dead_volume = 0.
                    if not hasattr(fld, 'tracer_object'):
                        fld.tracer_object = None
                    if fld.effluent.cp_normalized == 0:
                        fld.effluent.cp_normalized = False
                    else:
                        fld.effluent.cp_normalized = True
                    if fld.effluent.abs_normalized == 0:
                        fld.effluent.abs_normalized = False
                    else:
                        fld.effluent.abs_normalized = True
                for mf in expt.multi_floods_list.signal_lists[0].objects:
                    vols = mf.effluent.volumes
                    s = vols.shape
                    if s[0] == 3:
                        new_vols = np.zeros((4, s[1]))
                        new_vols[1, :] = np.nan
                        new_vols[0, :] = vols[0, :]
                        new_vols[2:, :] = vols[1:, :]
                        mf.effluent.volumes = new_vols
                    if not hasattr(mf.effluent, 'swelling_factor'):
                        mf.effluent.swelling_factor = 1.
                    if not hasattr(mf.effluent, 'revisions'):
                        mf.effluent.revisions = np.zeros((4, s[1]))
                        mf.effluent.revisions[1, :] = np.nan
                    if not hasattr(mf.effluent, 'micro_emulsions'):
                        mf.effluent.micro_emulsions = np.zeros((2, s[1]))
                        mf.effluent.micro_emulsions[:, :] = np.nan
                    if not hasattr(mf.effluent, 'ions') or not mf.effluent.ions.size:
                        mf.effluent.ions = np.zeros((10, s[1]))
                        mf.effluent.ions[:, :] = np.nan
                    if not hasattr(mf.effluent, 'dead_volume'):
                        mf.effluent.dead_volume = 0.
                    if not hasattr(mf, 'tracer_object'):
                        mf.tracer_object = None
                    if mf.effluent.cp_normalized == 0:
                        mf.effluent.cp_normalized = False
                    else:
                        mf.effluent.cp_normalized = True
                    if mf.effluent.abs_normalized == 0:
                        mf.effluent.abs_normalized = False
                    else:
                        mf.effluent.abs_normalized = True

        if version_list[i] == '0.6.0':
            for expt in self.flood_experiment_list.signal_lists[0].objects:

                for fld in expt.floods:
                    if not hasattr(fld.effluent, 'cp_tracer_object'):
                        fld.effluent.cp_tracer_object = None
                    if not hasattr(fld.effluent, 'abs_tracer_object'):
                        fld.effluent.abs_tracer_object = None
                    if fld.tracer_object is not None:
                        if not hasattr(fld.tracer_object, 'non_nan'):
                            fld.tracer_object = None
                    if fld.tracer_object is not None and not hasattr(fld.tracer_object, '_ipv'):
                        fld.tracer_object.__dict__['_ipv'] = 0.

                for mf in expt.multi_floods_list.signal_lists[0].objects:
                    if not hasattr(mf.effluent, 'cp_tracer_object'):
                        mf.effluent.cp_tracer_object = None
                    if not hasattr(mf.effluent, 'abs_tracer_object'):
                        mf.effluent.abs_tracer_object = None
                    if mf.tracer_object is not None:
                        if not hasattr(mf.tracer_object, 'non_nan'):
                            mf.tracer_object = None

        if version_list[i] == '0.6.1':
            for expt in self.flood_experiment_list.signal_lists[0].objects:

                if not hasattr(expt, 'notes'):
                    expt.notes = ''
                if not hasattr(expt, 'objective'):
                    expt.objective = ''

                for fld in expt.floods:
                    if not hasattr(fld, 'notes'):
                        fld.notes = ''
                    if not hasattr(fld, 'observations'):
                        fld.observations = ''
                    if not hasattr(fld.effluent, 'retention_ppm_min'):
                        fld.effluent.retention_ppm_min = 0.

                for mf in expt.multi_floods_list.signal_lists[0].objects:
                    if not hasattr(mf, 'notes'):
                        mf.notes = ''
                    if not hasattr(mf, 'observations'):
                        mf.observations = ''
                    if not hasattr(mf.effluent, 'retention_ppm_min'):
                        mf.effluent.retention_ppm_min = 0.

        if version_list[i] == '0.6.2':
            pass

        if version_list[i] == '0.6.3':
            pass

        if version_list[i] == '0.7.0':
            pass

        if version_list[i] == '0.7.1':

            for expt in self.flood_experiment_list.signal_lists[0].objects:
                for inj_fld in expt.injection_fluids_list.signal_lists[0].objects:
                    if isinstance(inj_fld.specific_fluid, fluid.BrineInjectionFluid):
                        if not hasattr(inj_fld.specific_fluid, 'stocks_list'):
                            inj_fld.specific_fluid.initialize_stocks_list()

        if version_list[i] == '0.8.0':

            for expt in self.flood_experiment_list.signal_lists[0].objects:

                for fld in expt.floods:
                    if not hasattr(fld, 'kro'):
                        fld.kro = np.nan

                for mf in expt.multi_floods_list.signal_lists[0].objects:
                    if not hasattr(mf, 'kro'):
                        mf.kro = np.nan

        if version_list[i] == '0.8.1':

            for oil in self.oil_list.signal_lists[0].objects:
                if not hasattr(oil, 'view_class'):
                    oil.view_class = fluid.OilTool
                if not hasattr(oil, 'api_gravity'):
                    oil.api_gravity = np.nan
                if not hasattr(oil, 'rho_m'):
                    oil.rho_m = np.nan
                if not hasattr(oil, 't_m'):
                    oil.t_m = np.nan
                if not hasattr(oil, 'gor'):
                    oil.gor = np.nan
                if not hasattr(oil, 'bubble_point'):
                    oil.bubble_point = np.nan
                if not hasattr(oil, 'gas_gravity'):
                    oil.gas_gravity = np.nan
                if not hasattr(oil, 'gas_mix'):
                    oil.gas_mix = fluid.GasMix([], [])
                if not hasattr(oil, 'rho_vs_p_data'):
                    oil.rho_vs_p_data = np.array([[], []])
                if not hasattr(oil, 'slm'):
                    oil.slm = utils.SignalListManager()
                    oil.slm.add_signal_list()
                    oil.slm.signal_lists[0].set_ref_lists([self.gas_list.signal_lists[0]])
                    oil.slm.signal_lists[0].objects = [oil.gas_mix]

            for oil_sample in self.oil_list.signal_lists[1].objects:
                if not hasattr(oil_sample, 'view_class'):
                    oil_sample.view_class = fluid.OilSampleTool
                if not hasattr(oil_sample, 'rheologies_list'):
                    oil_sample.rheologies_list = utils.SignalListManager()
                    oil_sample.rheologies_list.add_signal_list()
                if not hasattr(oil_sample, 'additives'):
                    oil_sample.additives = []
                if not hasattr(oil_sample, 'concentrations'):
                    oil_sample.concentrations = []
                if not hasattr(oil_sample, 'composition_set'):
                    oil_sample.composition_set = False

            for expt in self.flood_experiment_list.signal_lists[0].objects:
                for fld in expt.floods:
                    if not hasattr(fld.effluent, 'read_temperature'):
                        fld.effluent.read_temperature = fld.temperature
                    if not hasattr(fld.effluent, 'flood'):
                        fld.effluent.flood = fld
                for mf in expt.multi_floods_list.signal_lists[0].objects:
                    if not hasattr(mf.effluent, 'read_temperature'):
                        mf.effluent.read_temperature = mf.temperature
                    if not hasattr(mf.effluent, 'flood'):
                        mf.effluent.flood = mf

        if version_list[i] == '0.9.0':
            pass

        if version_list[i] == '0.9.1':
            pass

        if version_list[i] == '0.9.2':
            pass

        if version_list[i] == '0.9.3':

            for brine in self.brine_list.signal_lists[0].objects:
                nacl = "NaCl"
                if nacl not in brine.additive_salts.keys():
                    brine.additive_salts[nacl] = [[1, 1], ["Na+", "Cl-"], 58.443]

                if nacl not in brine.add_salt_composition.keys():
                    brine.add_salt_composition[nacl] = 0.

        if version_list[i] == '0.9.4':
            pass

        if version_list[i] == '0.9.5':
            pass

        if version_list[i] == '0.9.6':
            pass

        if version_list[i] == '0.9.7':

            print('upgrading to version 0.10.0')

            for oil in self.oil_list.signal_lists[0].objects:
                if not hasattr(oil, 'ift'):
                    oil.ift = np.nan
            for oil_sample in self.oil_list.signal_lists[1].objects:
                if not hasattr(oil_sample, 'ift'):
                    oil_sample.ift = np.nan
            for expt in self.flood_experiment_list.signal_lists[0].objects:
                if not hasattr(expt, 'ift_dict'):
                    expt.ift_dict = dict()

    @staticmethod
    def load_list_widget(lw: utils.SignalListManagerWidget, p_name_list: list):

        for i, name_list in enumerate(p_name_list):
            for name in name_list:
                lw.lists[i].addItem(name)

    def closeEvent(self, a0: QCloseEvent):
        already = False
        if self.source_file is not None:
            already = osp.exists(self.source_file.name)

        if already:
            dlg = QDialog(parent=self)
            lyt = QGridLayout()
            dlg.setLayout(lyt)
            lyt.addWidget(QLabel(parent=dlg, text='Overwrite?'), 0, 0, 1, 2)
            yes = QPushButton(parent=dlg, text='Yes')
            yes.clicked.connect(self.dump_pickle)
            yes.clicked.connect(dlg.close)
            no = QPushButton(parent=dlg, text='No')
            no.clicked.connect(dlg.close)
            no.clicked.connect(self.release_lock)
            lyt.addWidget(yes, 1, 0, 1, 1)
            lyt.addWidget(no, 1, 1, 1, 1)
            dlg.exec()
        else:
            msg = QMessageBox(self, text='Dumping to ' + self.project_manager.name + '.cfm')
            msg.exec_()
            self.dump_pickle()

    def dump_pickle(self, from_save_button=False):

        pickle_obj = self.create_pickle()

        if self.source_file is not None:
            self.release_lock()
            f_name = self.source_file.name
            copyfile(f_name, f_name[:-4] + '_copy_.cfm')
            self.source_file.close()
            self.source_file = None
            try:
                pickle_out = open(f_name, 'wb')
                pkl.dump(pickle_obj, pickle_out)
                QMessageBox(parent=self, text=f_name + ' saved.').exec_()
                pickle_out.close()
            except Exception as e:
                print('exception!')
                QMessageBox(parent=self, text=str(e)).exec_()
                copyfile(f_name[:-4] + '_copy_.cfm', f_name)

            print('removing')
            try:

                remove(f_name[:-4] + '_copy_.cfm')
            except Exception as e:
                print(e)
                # QMessageBox(parent=self, text=str(e)).exec_()
            print('removed')
            if from_save_button:
                self.source_file = open(f_name, 'rb')
                self.file_lock = FileLock(f_name + '.lock')
                self.file_lock.acquire()
        else:
            pickle_out = open(self.project_manager.name + '.cfm', 'wb')
            pkl.dump(pickle_obj, pickle_out)
            pickle_out.close()

    def release_lock(self):
        self.file_lock.release()
        self.file_lock = None


class ProjectManagerSaveStructure:

    def __init__(self, pmv: ProjectManagerView):

        pml_attr = 'project_manager_list'

        self.project_manager = pmv.project_manager
        # self.version = pmv.project_manager.version
        # print(self.version)
        self.plug_list = pmv.plug_list
        for c_list in self.plug_list.signal_lists:
            for obj in c_list.objects:
                obj.project_manager_list = None
        self.plug_name_list = self.extract_names(pmv.plug_list_widget)

        self.core_list = pmv.core_list
        for c_list in self.core_list.signal_lists:
            for obj in c_list.objects:
                obj.project_manager_list = None
        self.core_name_list = self.extract_names(pmv.core_list_widget)

        self.oil_list = pmv.oil_list
        for c_list in self.oil_list.signal_lists:
            for obj in c_list.objects:
                obj.project_manager_list = None
                if c_list == self.oil_list.signal_lists[1]:
                    for rheo in obj.rheologies_list.signal_lists[0].objects:
                        rheo.project_manager_list = None
        self.oil_name_list = self.extract_names(pmv.oil_list_widget)

        self.flood_experiment_list = pmv.flood_experiment_list
        for c_list in self.flood_experiment_list.signal_lists:
            for obj in c_list.objects:
                obj.project_manager_list = None
                for if_obj in obj.injection_fluids_list.signal_lists[0].objects:
                    if hasattr(if_obj, pml_attr):
                        if_obj.project_manager_list = None
                for fld_obj in obj.floods_list.signal_lists[0].objects:
                    if hasattr(fld_obj, pml_attr):
                        fld_obj.project_manager_list = None
                for mcf_obj in obj.multi_floods_list.signal_lists[0].objects:
                    if hasattr(mcf_obj, pml_attr):
                        mcf_obj.project_manager_list = None
                    if hasattr(mcf_obj, 'flood_view'):
                        mcf_obj.flood_view = None
                for c_flood in obj.floods:
                    c_flood.flood_view = None
        self.flood_experiment_name_list = self.extract_names(pmv.flood_experiment_list_widget)

        self.brine_list = pmv.brine_list
        for c_list in self.brine_list.signal_lists:
            for obj in c_list.objects:
                obj.project_manager_list = None
        self.brine_name_list = self.extract_names(pmv.brine_list_widget)

        self.chemical_list = pmv.chemical_list
        self.rheo_name_list = []
        for c_list in self.chemical_list.signal_lists:
            for obj in c_list.objects:
                if obj is not None:
                    obj.project_manager_list = None
                    if isinstance(obj, fluid.Polymer) and hasattr(obj, 'rheologies_list'):
                        # print(obj.__dict__)
                        for rheo in obj.rheologies_list.signal_lists[0].objects:
                            rheo.project_manager_list = None
                        #     rheo.base_fluid.project_manager_list = None
                    if isinstance(obj, fluid.Polymer) and hasattr(obj, 'base_fluids_list'):
                        # print(obj.__dict__)
                        for bf in obj.base_fluids_list.signal_lists[0].objects:
                            bf.project_manager_list = None
                        #     rheo.base_fluid.project_manager_list = None
        self.chemical_name_list = self.extract_names(pmv.chemical_list_widget)

        self.gas_list = pmv.gas_list
        for c_list in self.gas_list.signal_lists:
            for obj in c_list.objects:
                obj.project_manager_list = None
        self.gas_name_list = self.extract_names(pmv.gas_list_widget)

    @staticmethod
    def extract_names(list_widget: utils.SignalListManagerWidget):

        final_name_list = []

        for c_list in list_widget.lists:
            name_list = []
            count = c_list.count()
            for i in range(count):
                item = c_list.item(i)
                name_list.append(item.text())
            final_name_list.append(name_list)

        return final_name_list


if __name__ == "__main__":

    import sys

    app = QApplication(sys.argv)
    form = ProjectManagerView()
    # print(form.pos())

    sys.exit(app.exec_())
