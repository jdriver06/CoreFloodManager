
from sys import argv, exit as sys_exit
from PyQt5.QtWidgets import QApplication, QFileDialog, QWidget, QMessageBox, QDialog, QPushButton, \
    QVBoxLayout, QHBoxLayout, QLineEdit, QLabel, QTableWidget, QTableWidgetItem, QComboBox, QFrame
from PyQt5.QtWidgets import QMainWindow, QCheckBox, QListWidget, QListWidgetItem
from PyQt5.QtChart import QChart, QChartView, QLineSeries, QLogValueAxis, QScatterSeries
from PyQt5.Qt import QGridLayout, QIcon, QGroupBox, QPointF, pyqtSignal
from PyQt5.QtCore import Qt, QObject
from PyQt5.QtGui import QKeyEvent
from PyQt5 import QtGui
from numpy import average, array, zeros, nan, isnan, where, unique, power, divide, exp, sqrt, log10, linspace, \
    floor, ceil, min as np_min, max as np_max, meshgrid, ones, multiply, logspace, round as np_round, ndarray, delete, \
    float as np_float, interp, sum as np_sum, logical_and, argsort
from pandas import read_excel, read_csv
from enum import Enum, auto
from scipy.optimize import curve_fit
from scipy.stats import linregress
from matplotlib import pyplot as plt, cm
import j_utils as utils
from math import log
from copy import copy, deepcopy
import reports


__version__ = '0.1.0'


class Fluid:

    def __init__(self):

        self.file = ""
        self.dir = ""
        self.name = ""

    def set_file(self):

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(None, "QFileDialog.getOpenFileName()", self.dir,
                                                         "CSV Files (*.csv);; EXCEL Files (*.xlsx)", options=options)
        if file_name:
            self.file = file_name
            i = file_name.rfind('/')
            self.dir = file_name[:i+1]


class FluidView(QDialog):

    def __init__(self, parent, *args):
        super(FluidView, self).__init__(parent=parent)

        sf = parent.parent().width() / 280.
        self.setFixedSize(int(sf * 200), int(sf * 50))
        self.setWindowTitle(' ')
        layout = QHBoxLayout()
        self.setLayout(layout)

        self.brine_button = QPushButton(parent=self, text='Brine-based')
        self.brine_button.clicked.connect(self.brine_button_clicked)
        self.brine_button.clicked.connect(self.close)
        self.oil_button = QPushButton(parent=self, text='Oil-based')
        self.oil_button.clicked.connect(self.oil_button_clicked)
        self.oil_button.clicked.connect(self.close)

        layout.addWidget(self.brine_button)
        layout.addWidget(self.oil_button)

        if isinstance(parent.parent(), PolymerTool):
            self.brine_button_clicked()
            self.close()
        else:
            self.show()

    def brine_button_clicked(self):
        try:
            BrineInjectionFluidView(self.parent())
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()

    def oil_button_clicked(self):
        FluidViewOil(self.parent())


class FluidViewOil(QDialog):

    def __init__(self, parent, *args):
        super(FluidViewOil, self).__init__(parent=parent)

        sf = parent.parent().width() / 280.
        self.setFixedSize(int(sf * 200), int(sf * 50))
        self.setWindowTitle(' ')
        layout = QHBoxLayout()
        self.setLayout(layout)

        self.dead_button = QPushButton(parent=self, text='Dead oil')
        self.dead_button.clicked.connect(self.dead_button_clicked)
        self.dead_button.clicked.connect(self.close)
        self.live_button = QPushButton(parent=self, text='Live oil')
        self.live_button.clicked.connect(self.live_button_clicked)
        self.live_button.clicked.connect(self.close)

        layout.addWidget(self.dead_button)
        layout.addWidget(self.live_button)

        if isinstance(parent.parent(), PolymerTool):
            self.dead_button_clicked()
            self.close()
        else:
            self.show()

    def dead_button_clicked(self):
        try:
            OilInjectionFluidView(self.parent())
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()

    def live_button_clicked(self):
        view = OilInjectionFluidView(self.parent())
        view.fluid_type = FluidType.LIVE_OIL


class FluidType(Enum):

    OIL = auto()
    LIVE_OIL = auto()
    BRINE = auto()
    POLYMER_SOLUTION = auto()
    UNRECOGNIZED = auto()


class BaseFluidType(Enum):

    BRINE_BASED = auto()
    OIL_BASED = auto()
    UNRECOGNIZED = auto()


class InjectionFluid:

    def __init__(self, name: str, ref_objects: list, *args):
        # print('made it to InjectionFluid constructor.')
        print('InjectionFluid __init__', name, ref_objects, *args)

        self.view_class = InjectionFluidView
        # print(name, ref_objects, args)

        self.name = name
        self.project_manager_list = None
        self.experiment = None

        self.ref_objects = ref_objects
        args = args[0]
        # print('Injection Fluid args', args)
        self.concentrations = args[0]

        if args[1] == 'brine-based':
            self.specific_fluid = BrineInjectionFluid(self)

        elif args[1] == 'oil-based':
            self.specific_fluid = OilInjectionFluid(self)

    def copy_me(self) -> object:

        args = [self.concentrations]
        if self.get_base_fluid_type() == BaseFluidType.BRINE_BASED:
            args.append('brine-based')
        else:
            args.append('oil-based')
        args = [args]

        new_injection_fluid = InjectionFluid(self.name + ' (copy)', self.ref_objects, *args)

        return new_injection_fluid

    def get_viscosity(self, temp: float, shear: float) -> float:
        # print('getting viscosity in InjectionFluid')
        return self.specific_fluid.get_viscosity(temp, shear)

    def get_density(self, temp: float) -> float:
        return self.specific_fluid.get_density(temp)

    def get_fluid_type(self) -> FluidType:

        if isinstance(self.specific_fluid, BrineInjectionFluid):
            if self.specific_fluid.polymer_solution is None:
                return FluidType.BRINE
            else:
                return FluidType.POLYMER_SOLUTION

        elif isinstance(self.specific_fluid, OilInjectionFluid):
            if self.specific_fluid.polymer_solution is None:
                return FluidType.LIVE_OIL
            else:
                return FluidType.OIL

        else:
            return FluidType.UNRECOGNIZED

    def get_base_fluid_type(self) -> BaseFluidType:

        if isinstance(self.specific_fluid, BrineInjectionFluid):
            return BaseFluidType.BRINE_BASED

        elif isinstance(self.specific_fluid, OilInjectionFluid):
            return BaseFluidType.OIL_BASED

        else:
            return BaseFluidType.UNRECOGNIZED


class InjectionFluidView(QObject, utils.SignalView):

    def __init__(self, injection_fluid: InjectionFluid, parent):
        super(InjectionFluidView, self).__init__()

        try:
            if hasattr(parent.list_widget.parent().parent(), 'experiment'):
                injection_fluid.experiment = parent.list_widget.parent().parent().experiment
        except Exception as e:
            QMessageBox(parent=parent, text=str(e)).exec_()

        if isinstance(injection_fluid.specific_fluid, BrineInjectionFluid):
            view = BrineInjectionFluidView(parent=parent.list_widget, brine_if=injection_fluid)
            view.close_signal.connect(self.relay_close_signal)
        elif isinstance(injection_fluid.specific_fluid, OilInjectionFluid):
            view = OilInjectionFluidView(parent=parent.list_widget, oil_if=injection_fluid)
            view.close_signal.connect(self.relay_close_signal)
            # if injection_fluid.specific_fluid.polymer_solution is not None:
            #     PolymerSolutionView(injection_fluid.specific_fluid.polymer_solution, parent=parent.list_widget)

    def relay_close_signal(self, _):
        self.close_signal.emit(self)


class BrineInjectionFluid:

    def __init__(self, injection_fluid: InjectionFluid):

        self.injection_fluid = injection_fluid
        self.brine = injection_fluid.ref_objects[0]
        self.stocks_list = {}
        self.stocks_base_fluids = {}
        self.mixing_mass = 500.
        self.formulation_mixing_masses = {}
        self.mixing_dilution_x_factor = 1.
        self.mixing_base_fluid = None
        self.mixing_dilution_fluid = None

        self.pH = nan
        self.RI = nan

        p_concentration = 0
        p = None

        if len(injection_fluid.ref_objects) > 1:
            non_brines = injection_fluid.ref_objects[1:]
        else:
            non_brines = []

        self.additives = non_brines

        for obj, concentration in zip(non_brines, injection_fluid.concentrations):
            if isinstance(obj, Polymer):
                p = obj
                p_concentration = concentration
                break

        if p_concentration > 0:
            self.polymer_solution = PolymerSolution(base_fluid=self.brine, primary_additive=p,
                                                    concentration=p_concentration)
        else:
            self.polymer_solution = None

        self.initialize_stocks_list()

    def initialize_stocks_list(self):

        self.stocks_list = {}
        self.stocks_base_fluids = {}
        self.mixing_mass = 500.
        self.formulation_mixing_masses = {}
        self.mixing_dilution_x_factor = 1.
        self.mixing_base_fluid = None
        self.mixing_dilution_fluid = None

        if self.brine.has_additive_salt():
            salts, concentrations = self.brine.additive_salts_and_concentrations()
            for s, c in zip(salts, concentrations):
                if c <= 10.:
                    self.stocks_list[s] = 10.
                else:
                    self.stocks_list[s] = c

        for additive in self.additives:
            if isinstance(additive, Surfactant) or isinstance(additive, CoSolvent):
                self.stocks_list[additive] = 100.
            if isinstance(additive, Formulation):
                self.stocks_list[additive] = 8.
                self.formulation_mixing_masses[additive] = 100.
                for item in additive.sc_list:
                    self.stocks_list[item] = 100.
            if isinstance(additive, Polymer):
                if self.polymer_solution.concentration <= 10000.:
                    self.stocks_list[additive] = 10000.
                else:
                    self.stocks_list[additive] = self.polymer_solution.concentration

    def get_viscosity(self, temp: float, shear=0.):

        ps = self.polymer_solution
        if ps is not None:
            if ps.stored_averages:
                return ps.get_viscosity(temp, shear)
            else:
                v = ps.primary_additive.estimate_viscosity(ps.base_fluid, temp, shear, ps.concentration)
                if not isnan(v):
                    return v

        return self.brine.get_viscosity(temp)

    def get_density(self, temp: float):
        return self.brine.get_density(temp)

    def has_additive_type(self, cls) -> bool:

        for additive in self.additives:
            if isinstance(additive, cls):
                return True

        return False

    def n_additives_of_type(self, cls) -> int:

        if not self.has_additive_type(cls):
            return 0

        n_adds_of_type = 0
        for additive in self.additives:
            if isinstance(additive, cls):
                n_adds_of_type += 1

        return n_adds_of_type

    def additives_and_concentrations(self, cls):

        if not self.has_additive_type(cls):
            return [], []

        adds = []
        concentrations = []

        for add, c in zip(self.additives, self.injection_fluid.concentrations):
            if isinstance(add, cls):
                adds.append(add)
                concentrations.append(c)

        return adds, concentrations

    def has_formulation(self) -> bool:

        return self.has_additive_type(Formulation)

    def n_formulations(self) -> int:

        return self.n_additives_of_type(Formulation)

    def formulations_and_concentrations(self):

        return self.additives_and_concentrations(Formulation)

    def has_surfactant(self) -> bool:

        return self.has_additive_type(Surfactant)

    def surfactants_and_concentrations(self):

        return self.additives_and_concentrations(Surfactant)

    def has_cosolvent(self) -> bool:

        return self.has_additive_type(CoSolvent)

    def cosolvents_and_concentrations(self):

        return self.additives_and_concentrations(CoSolvent)

    def has_polymer(self) -> bool:

        if self.polymer_solution is not None:
            return True

        return False

    def n_polymers(self) -> int:

        return self.n_additives_of_type(Polymer)

    def polymers_and_concentrations(self):

        return self.additives_and_concentrations(Polymer)

    def has_alkali(self) -> bool:

        return self.brine.has_alkali()

    def has_additive_salt(self) -> bool:

        return self.brine.has_additive_salt()

    def has_additive_non_salt(self) -> bool:

        return self.brine.has_additive_non_salt()

    def scavengers_and_concentrations(self):

        return self.additives_and_concentrations(Scavenger)

    def others_and_concentrations(self):

        return self.additives_and_concentrations(Other)


class Oil:

    def __init__(self, name: str, pm=None):

        self.name = name
        self.project_manager_list = pm
        self.view_class = OilTool

        self.api_gravity = nan
        self.ift = nan
        self.rho_m = nan
        self.t_m = nan
        self.gor = nan
        self.bubble_point = nan
        self.gas_gravity = nan
        self.gas_mix = GasMix([], [])
        self.slm = utils.SignalListManager()
        self.slm.add_signal_list()
        self.slm.signal_lists[0].objects = [self.gas_mix]
        self.rho_vs_p_data = array([[], []])

    def get_rho_60(self) -> float:

        if isnan(self.api_gravity):
            return nan

        sg = 141.5 / (131.5 + self.api_gravity)
        return 0.99907 * sg

    def get_alpha_60(self) -> float:

        if isnan(self.api_gravity) or isnan(self.rho_m) or isnan(self.t_m):
            return nan

        return (sqrt(1. + 3.2 * log(self.get_rho_60() / self.rho_m)) - 1.) / (1.6 * (self.t_m - 60.))

    def get_rho_v_t(self, t: array) -> array:

        alpha60 = self.get_alpha_60()
        dt = t - 60.
        cTL = exp(multiply(-alpha60 * dt, 1. + 0.8 * alpha60 * dt))

        return self.get_rho_60() * cTL

    def get_density(self, t: float) -> float:

        t = t * (9. / 5.) + 32.
        alpha60 = self.get_alpha_60()
        dt = t - 60.
        cTL = exp(multiply(-alpha60 * dt, 1. + 0.8 * alpha60 * dt))

        return self.get_rho_60() * cTL

    def get_live_oil_density(self, p: float) -> float:

        p_data = self.rho_vs_p_data[0, :]
        p_data.flatten()
        rho_data = self.rho_vs_p_data[1, :]
        rho_data.flatten()

        if not p_data.size:
            return nan

        elif p_data.size == 1:
            return rho_data[0]

        elif p < p_data[0]:
            return rho_data[0]

        elif p > p_data[-1]:
            return rho_data[-1]

        return interp(p, p_data, rho_data)

    def get_live_oil_composition(self):

        # gas density [kg/m3]
        rho_g = self.gas_gravity * 1.222
        # gas / oil volume ratio [m3/m3]
        volume_ratio = self.gor * 0.17810760667903525
        # oil density [kg/m3]
        rho_o = self.get_rho_60() * 1.e3

        # gas / oil mass ratio [kg/kg]
        mass_ratio = volume_ratio * (rho_g / rho_o)

        ppm = []
        for c in self.gas_mix.concentrations:
            ppm.append(mass_ratio * c * 1.e4)

        print('gas/oil mass ratio: {:.6f}'.format(mass_ratio))

        return self.gas_mix.gases, ppm


class OilTool(QDialog, utils.SignalView):

    editValueChanged = pyqtSignal()

    def __init__(self, oil: Oil, parent):
        super(OilTool, self).__init__(parent=parent)

        sf = parent.sf
        self.sf = sf
        self.oil = oil
        if not oil.slm.signal_lists[0].ref_lists:
            oil.slm.signal_lists[0].set_ref_lists([parent.gas_list.signal_lists[0]])

        self.setWindowTitle('{}'.format(oil.name))
        self.setFixedSize(int(450 * sf), int(190 * sf))
        lyt = QGridLayout()
        self.setLayout(lyt)

        api_label = QLabel(parent=self, text=u'API gravity [' + chr(176) + u']:')
        api_edit = QLineEdit(parent=self)
        api_edit.editingFinished.connect(self.api_edited)
        ift_label = QLabel(parent=self, text='IFT [dynes/cm]:')
        ift_edit = QLineEdit(parent=self)
        ift_edit.editingFinished.connect(self.ift_edited)
        m_label = QLabel(parent=self, text='Measured Density at Temperature')
        rho_label = QLabel(parent=self, text=chr(961) + u' [g/cc]:')
        rho_edit = QLineEdit(parent=self)
        rho_edit.editingFinished.connect(self.rho_edited)
        t_label = QLabel(parent=self, text=u'T [' + chr(176) + u'F]:')
        t_edit = QLineEdit(parent=self)
        t_edit.editingFinished.connect(self.t_edited)
        alpha_label = QLabel(parent=self, text='')

        gor_label = QLabel(parent=self, text='GOR [scf/bbl]:')
        gor_edit = QLineEdit(parent=self)
        gor_edit.editingFinished.connect(self.gor_edited)
        bp_label = QLabel(parent=self, text='Bubble Point [psi]:')
        bp_edit = QLineEdit(parent=self)
        bp_edit.editingFinished.connect(self.bp_edited)
        gravity_label = QLabel(parent=self, text='Gas Gravity:')
        gravity_edit = QLineEdit(parent=self)
        gravity_edit.editingFinished.connect(self.gravity_edited)

        self.show_density_button = QPushButton(parent=self, text=u'Show ' + chr(961) + u' vs. T')
        self.show_density_button.setEnabled(False)
        self.show_density_button.clicked.connect(self.show_density)

        self.gas_composition_button = QPushButton(parent=self, text=' Gas Composition \n[%wt]')
        self.gas_composition_button.clicked.connect(self.view_edit_gas_composition)
        self.gas_composition_text = QLabel(parent=self, text='')
        font = self.gas_composition_text.font()
        font.setFamily('Courier')
        self.gas_composition_text.setFont(font)
        self.gas_composition_text.setAlignment(Qt.AlignTop)
        self.gas_composition_text.setStyleSheet('padding :1px')
        self.live_density_button = QPushButton(parent=self, text=u'L.O. ' + chr(961) + u' vs. P')
        self.live_density_button.clicked.connect(self.view_edit_rho_vs_p_data)

        lyt.addWidget(api_label, 0, 0, 1, 1)
        lyt.addWidget(api_edit, 0, 1, 1, 1)
        lyt.addWidget(ift_label, 1, 0, 1, 1)
        lyt.addWidget(ift_edit, 1, 1, 1, 1)
        lyt.addWidget(m_label, 2, 0, 1, 2)
        lyt.addWidget(rho_label, 3, 0, 1, 1)
        lyt.addWidget(rho_edit, 3, 1, 1, 1)
        lyt.addWidget(t_label, 4, 0, 1, 1)
        lyt.addWidget(t_edit, 4, 1, 1, 1)
        lyt.addWidget(alpha_label, 5, 0, 1, 2)
        lyt.addWidget(self.show_density_button, 6, 0, 1, 2)
        lyt.addWidget(gor_label, 0, 2, 1, 1)
        lyt.addWidget(gor_edit, 0, 3, 1, 1)
        lyt.addWidget(bp_label, 1, 2, 1, 1)
        lyt.addWidget(bp_edit, 1, 3, 1, 1)
        lyt.addWidget(gravity_label, 2, 2, 1, 1)
        lyt.addWidget(gravity_edit, 2, 3, 1, 1)
        lyt.addWidget(self.gas_composition_button, 3, 2, 2, 1)
        lyt.addWidget(self.gas_composition_text, 3, 3, 2, 1)
        lyt.addWidget(self.live_density_button, 6, 2, 1, 2)

        self.api_edit = api_edit
        self.sync_edit(api_edit, oil.api_gravity)
        self.ift_edit = ift_edit
        self.sync_edit(ift_edit, oil.ift)
        self.rho_edit = rho_edit
        self.sync_edit(rho_edit, oil.rho_m)
        self.t_edit = t_edit
        self.sync_edit(t_edit, oil.t_m)
        self.gor_edit = gor_edit
        self.sync_edit(gor_edit, oil.gor)
        self.bp_edit = bp_edit
        self.sync_edit(bp_edit, oil.bubble_point)
        self.alpha_label = alpha_label
        self.gravity_edit = gravity_edit
        self.sync_edit(self.gravity_edit, oil.gas_gravity)

        self.editValueChanged.connect(self.update_alpha_60)
        self.editValueChanged.emit()
        self.update_gc_text()

        self.gas_mix_tool = None
        self.lo_density_tool = None

        self.show()

    @staticmethod
    def sync_edit(edit: QLineEdit, val: float):

        edit.blockSignals(True)

        if isnan(val):
            edit.setText('')
        else:
            edit.setText('{}'.format(val))

        edit.blockSignals(False)

    def api_edited(self):

        txt = self.api_edit.text()

        try:
            val = float(txt)
            self.oil.api_gravity = val
            self.sync_edit(self.api_edit, val)
            self.editValueChanged.emit()

        except ValueError:
            self.sync_edit(self.api_edit, self.oil.api_gravity)

    def ift_edited(self):

        txt = self.ift_edit.text()

        try:
            val = float(txt)
            if val <= 0.:
                raise ValueError

            self.oil.ift = val
            self.sync_edit(self.ift_edit, val)
            self.editValueChanged.emit()

        except ValueError:
            self.sync_edit(self.ift_edit, self.oil.ift)

    def rho_edited(self):

        txt = self.rho_edit.text()

        try:
            val = float(txt)
            if val <= 0.:
                raise ValueError

            self.oil.rho_m = val
            self.sync_edit(self.rho_edit, val)
            self.editValueChanged.emit()

        except ValueError:
            self.sync_edit(self.rho_edit, self.oil.rho_m)

    def t_edited(self):

        txt = self.t_edit.text()

        try:
            val = float(txt)
            self.oil.t_m = val
            self.sync_edit(self.t_edit, val)
            self.editValueChanged.emit()

        except ValueError:
            self.sync_edit(self.t_edit, self.oil.t_m)

    def gor_edited(self):

        txt = self.gor_edit.text()

        try:
            val = float(txt)
            self.oil.gor = val
            self.sync_edit(self.gor_edit, val)
            self.editValueChanged.emit()

        except ValueError:
            self.sync_edit(self.gor_edit, self.oil.gor)

    def bp_edited(self):

        txt = self.bp_edit.text()

        try:
            val = float(txt)
            self.oil.bubble_point = val
            self.sync_edit(self.bp_edit, val)
            self.editValueChanged.emit()

        except ValueError:
            self.sync_edit(self.bp_edit, self.oil.bubble_point)

    def gravity_edited(self):

        txt = self.gravity_edit.text()

        try:
            val = float(txt)
            self.oil.gas_gravity = val
            self.sync_edit(self.gravity_edit, val)
            self.editValueChanged.emit()

        except ValueError:
            self.sync_edit(self.gravity_edit, self.oil.gas_gravity)

    def update_alpha_60(self):

        alpha60 = self.oil.get_alpha_60()
        self.alpha_label.setText(chr(945) + u'60 = {:.6f} [cc/g-F]'.format(alpha60))

        if not isnan(alpha60):
            self.show_density_button.setEnabled(True)
        else:
            self.show_density_button.setEnabled(False)

    def show_density(self):

        t = array([i for i in range(60, 213)])
        rho_v_t = self.oil.get_rho_v_t(t)

        plt.figure()
        plt.plot((t - 32.) * (5. / 9.), rho_v_t)
        plt.xlabel(u'T [' + chr(176) + u'C]')
        plt.ylabel(chr(961) + u' [g/cc]')
        plt.title('{} Density vs. Temperature'.format(self.oil.name))
        plt.show()

    def view_edit_gas_composition(self):
        if self.gas_mix_tool is None:
            self.gas_mix_tool = GasMixView(self.oil.gas_mix, self)

    def view_edit_rho_vs_p_data(self):
        if self.lo_density_tool is None:
            self.lo_density_tool = LiveOilDensityVsTDialog(self)

    def update_gc_text(self):

        max_len = 0

        for gas in self.oil.gas_mix.gases:
            name_len = len(gas.name)
            if name_len > max_len:
                max_len = name_len

        gc_text = ''

        i = 1
        for gas, c in zip(self.oil.gas_mix.gases, self.oil.gas_mix.concentrations):
            if i < 3:
                name_text = '{:>' + str(max_len) + '}: '
                gc_text += name_text.format(gas.name)
                gc_text += '{:5.1f}\n'.format(c)
            else:
                name_text = '{:<' + str(max_len) + '}'
                gc_text += name_text.format('...')
                break
            i += 1

        self.gas_composition_text.setText(gc_text)
        self.oil.get_live_oil_composition()

    def closeEvent(self, a0: QtGui.QCloseEvent):

        if self.gas_mix_tool is not None:
            self.gas_mix_tool.close()
        if self.lo_density_tool is not None:
            self.lo_density_tool.close()

        self.close_signal.emit(self)


class LiveOilDensityVsTDialog(QDialog):

    def __init__(self, parent: OilTool):
        super(LiveOilDensityVsTDialog, self).__init__(parent=parent)
        self.setWindowTitle('L.O. Density')
        self.sf = parent.sf
        self.oil_tool = parent
        self.setFixedSize(int(self.sf * 250), int(self.sf * 300))
        self.setLayout(QVBoxLayout())
        self.table = utils.SignalTable(self, parent.oil.rho_vs_p_data)
        for c in range(self.table.columnCount()):
            self.table.thresholds[c] = 0.
        self.table.setHorizontalHeaderLabels(['P [psi]', chr(961) + u' [g/cc]'])
        self.layout().addWidget(self.table)

        self.show()

    def closeEvent(self, a0: QtGui.QCloseEvent):

        data = self.table.data
        p_data = data[0, :]
        rho_data = data[1, :]

        non_nans = logical_and(~isnan(p_data), ~isnan(rho_data))
        p_data = p_data[non_nans]
        p_data.flatten()
        rho_data = rho_data[non_nans]
        rho_data.flatten()

        if not p_data.size:
            self.oil_tool.oil.rho_vs_p_data = array([[], []])
        else:
            ind = argsort(p_data)
            p_data = p_data[ind]
            rho_data = rho_data[ind]
            self.oil_tool.oil.rho_vs_p_data = array([p_data.tolist(), rho_data.tolist()])

        self.oil_tool.lo_density_tool = None


class OilSample:

    def __init__(self, name: str, ref_oils: list):
        self.name = name
        self.ref_objects = ref_oils
        self.ift = nan
        self.project_manager_list = None
        self.view_class = OilSampleTool
        self.rheologies_list = utils.SignalListManager()
        self.rheologies_list.add_signal_list()

        self.additives = []
        self.concentrations = []
        self.composition_set = False

    def get_ift(self) -> float:

        if not isnan(self.ift):
            return self.ift

        if isinstance(self.ref_objects[0], Oil):
            return self.ref_objects[0].ift

        return nan


class Diluent:

    def __init__(self, name: str, pm=None):
        self.name = name
        self.project_manager_list = pm


class OilSampleTool(QDialog, utils.SignalView):

    def __init__(self, oil_sample: OilSample, parent, *args):
        super(OilSampleTool, self).__init__(parent=parent)
        self.setLayout(QHBoxLayout())
        self.oil_sample = oil_sample
        sf = parent.sf
        self.sf = sf
        self.setFixedSize(int(sf * 300), int(sf * 300))
        self.setWindowTitle(oil_sample.name + ' Viewer')

        list_name = ['Rheologies']
        cls = [[oil_sample_rheology_view_wrapper], [InjectionFluid]]
        args = [['Name'], [oil_sample]]
        self.slm = self.oil_sample.rheologies_list
        self.list_widget = utils.SignalListManagerWidget(self, self.oil_sample.rheologies_list, list_name, cls, *args)
        self.slm.add_signal_list()
        for rheo in oil_sample.rheologies_list.signal_lists[0].objects:
            rheo.project_manager_list = self.list_widget

        self.list_widget.pm_list.signal_lists[0].set_ref_lists([parent.oil_list.signal_lists[1],
                                                                parent.oil_list.signal_lists[2]])

        self.layout().addWidget(self.list_widget)

        self.show()
        self.composition_view = None
        self.ift_view = None
        self.view_composition()

        if not oil_sample.composition_set:
            self.list_widget.setEnabled(False)
        elif len(oil_sample.additives) > 0:
            self.view_ift()

    def view_composition(self):

        self.composition_view = OilInjectionFluidView(parent=self)
        self.composition_view.enable_combo_box = False
        self.composition_view.setWindowTitle('Sample Composition: ' + self.oil_sample.name)
        self.composition_view.name_text.setText('Base Oil: ' + self.oil_sample.ref_objects[0].name)
        self.composition_view.fluid_name_edit.setText(self.oil_sample.ref_objects[0].name)
        self.composition_view.fluid_name_edit.setEnabled(False)
        # self.composition_view.name_text.setVisible(False)
        self.composition_view.fluid_name_edit.setVisible(False)
        i = self.parent().oil_list.signal_lists[1].objects.index(self.oil_sample)
        self.composition_view.source_fluid_combo_box.setCurrentIndex(i)
        self.composition_view.source_fluid_combo_box.setEnabled(False)
        self.composition_view.source_fluid_combo_box.setVisible(False)
        self.composition_view.submit_button.clicked.connect(self.enable_rheology_list)

        self.composition_view.setFixedHeight(int(70 * self.sf))

        self.composition_view.add_button.setEnabled(True)
        for add, c in zip(self.oil_sample.additives, self.oil_sample.concentrations):
            print(add, c)
            self.composition_view.add_button.click()
            w = self.composition_view.combo_widgets[-1]
            w.set_combobox_to_name(add.name)
            w.combo_box_2.setEnabled(False)
            w.edit.setText(str(c))
            w.edit.setEnabled(False)
            print('additive added.')

        if self.oil_sample.composition_set:
            self.composition_view.add_button.setEnabled(False)
            self.composition_view.remove_button.setEnabled(False)
            self.composition_view.submit_button.setEnabled(False)
        else:
            self.composition_view.add_button.setEnabled(True)
            self.composition_view.remove_button.setEnabled(True)
            self.composition_view.submit_button.setEnabled(True)

    def view_ift(self):

        if self.ift_view is None:
            self.ift_view = OilSampleIFTTool(self)
            self.ift_view.show()

    def enable_rheology_list(self):

        self.list_widget.setEnabled(True)

    def closeEvent(self, a0: QtGui.QCloseEvent) -> None:
        super(OilSampleTool, self).closeEvent(a0)
        self.close_signal.emit(self)


class OilSampleIFTTool(QDialog):

    def __init__(self, parent: OilSampleTool):
        super(OilSampleIFTTool, self).__init__(parent=parent)
        self.parent = parent
        self.sf = parent.sf
        self.setWindowTitle('{} IFT Tool'.format(parent.oil_sample.name))
        self.setFixedSize(int(self.sf * 250), int(self.sf * 50))

        lyt = QHBoxLayout()
        self.setLayout(lyt)

        ift_label = QLabel(parent=self, text='IFT [dynes/cm]:')
        self.ift_edit = QLineEdit(parent=self)
        self.ift_edit.editingFinished.connect(self.ift_edited)

        lyt.addWidget(ift_label)
        lyt.addWidget(self.ift_edit)

        self.load_edit()

    def load_edit(self):

        self.ift_edit.setText('')
        ift = self.parent.oil_sample.ift
        if not isnan(ift):
            self.ift_edit.setText(str(ift))

    def ift_edited(self):

        txt = self.ift_edit.text()

        try:
            val = float(txt)
            if val < 0.:
                raise ValueError

            self.parent.oil_sample.ift = val

        except ValueError:
            self.load_edit()

    def closeEvent(self, a0: QtGui.QCloseEvent) -> None:

        try:
            self.parent.ift_view = None
        except Exception as e:
            print(e)


def oil_sample_rheology_view_wrapper(parent, *args):

    try:
        view = OilInjectionFluidView(parent=parent)
        oil_sample = args[1][0]

        i = parent.parent().parent().oil_list.signal_lists[1].objects.index(oil_sample)
        view.source_fluid_combo_box.setCurrentIndex(i)
        view.name_text.setText('Rheology Name:')
        view.setWindowTitle(oil_sample.name + ': Oil-based Rheology')
        view.enable_combo_box = False
        view.source_fluid_combo_box.setEnabled(False)

        if oil_sample.composition_set:
            # for add, c in zip(oil_sample.additives, oil_sample.concentrations):
            #     print(add, c)
            #     view.add_button.click()
            #     w = view.combo_widgets[-1]
            #     w.set_combobox_to_name(add.name)
            #     w.combo_box_2.setEnabled(False)
            #     w.edit.setText(str(c))
            #     w.edit.setEnabled(False)
            #     w.setVisible(False)
            view.add_button.setVisible(False)
            view.remove_button.setVisible(False)
            view.source_fluid_combo_box.setVisible(False)
            view.setFixedHeight(int(100 * parent.parent().sf))

    except Exception as e:
        QMessageBox(parent=parent, text=str(e)).exec_()


class Gas:

    def __init__(self, name: str, pm=None):
        self.name = name
        self.project_manager_list = pm


class GasMix:

    def __init__(self, gases: list, concentrations: list):
        self.ref_objects = gases
        self.concentrations = concentrations

    @property
    def gases(self):
        return self.ref_objects


class GasMixView(QDialog):

    def __init__(self, gas_mix: GasMix, parent):
        super(GasMixView, self).__init__(parent=parent)

        self.setFixedSize(int(parent.sf * 350), int(parent.sf * 70))
        self.setWindowTitle(parent.oil.name + ' Gas Mix')
        self.setLayout(QVBoxLayout())
        self.gas_mix = gas_mix
        self.sf = parent.sf

        button_widget = QWidget()
        button_widget.setLayout(QHBoxLayout())
        self.add_button = QPushButton(parent=self, text='+')
        self.add_button.setFixedSize(int(parent.sf * 25), int(parent.sf * 25))
        self.add_button.clicked.connect(self.add_gas)
        self.remove_button = QPushButton(parent=self, text='-')
        self.remove_button.setFixedSize(int(parent.sf * 25), int(parent.sf * 25))
        self.remove_button.clicked.connect(self.remove_gas)
        button_widget.layout().addWidget(self.add_button)
        button_widget.layout().addWidget(self.remove_button)
        self.layout().addWidget(button_widget)

        self.combo_widgets = []

        self.build_from_stored()
        self.show()

    def build_from_stored(self):

        for gas, c in zip(self.gas_mix.gases, self.gas_mix.concentrations):
            self.add_gas()
            widget = self.combo_widgets[-1]
            widget.set_combobox_to_name(gas.name)
            widget.edit.setText(str(c))

    def capture_from_stored(self):

        pmv = self.parent().parent()
        gases = []
        concentrations = []

        for widget in self.combo_widgets:
            txt = widget.edit.text()
            gas_name = widget.combo_box_2.currentText()
            if txt and gas_name:
                c = float(txt)
                concentrations.append(c)
                gas = pmv.gas_list.signal_lists[0].find_item_by_name(gas_name)
                gases.append(gas)

        self.gas_mix.ref_objects = gases
        if concentrations:
            concentrations = array(concentrations)
            concentrations = 100. * concentrations / np_sum(concentrations)
            self.gas_mix.concentrations = list(np_round(concentrations, 1))
        else:
            self.gas_mix.concentrations = []

    def add_gas(self, conc_index: int=-1):

        pmv = self.parent().parent()
        new_widget = AdditiveWidgetGasMix([pmv.oil_list, pmv.gas_list], self, conc_index)
        self.layout().addWidget(new_widget)
        self.combo_widgets.append(new_widget)
        self.setFixedHeight(self.height() + int(self.sf * 51))

    def remove_gas(self):

        if not self.combo_widgets:
            return

        widget = self.combo_widgets.pop()
        self.layout().removeWidget(widget)
        widget.destroy()
        self.setFixedHeight(self.height() - int(self.sf * 51))

    def closeEvent(self, _):
        self.capture_from_stored()
        self.parent().update_gc_text()
        self.parent().gas_mix_tool = None


class OilInjectionFluid:

    def __init__(self, injection_fluid: InjectionFluid):

        self.injection_fluid = injection_fluid
        self.oil_sample = injection_fluid.ref_objects[0]
        self._inferred_viscosity = nan

        self.additives = []
        if len(injection_fluid.ref_objects) > 1:
            self.additives = injection_fluid.ref_objects[1:]

        d_concentration = 0
        d = None

        non_oils = self.additives

        # print('OilInjectionFluid: ', non_oils, injection_fluid.concentrations)

        is_live_oil = False
        if non_oils:
            for obj, concentration in zip(non_oils, injection_fluid.concentrations):
                if isinstance(obj, Diluent) or isinstance(obj, Gas):
                    d = obj
                    d_concentration = concentration
                    if isinstance(obj, Gas):
                        is_live_oil = True
                    break

        # print('constructing PolymerSolution.')
        if not is_live_oil:
            self.polymer_solution = PolymerSolution(base_fluid=self.oil_sample, primary_additive=d,
                                                    concentration=d_concentration)
        else:
            self.polymer_solution = None

    def get_viscosity(self, temp: float, shear: float):
        if self.polymer_solution is not None:
            return self.polymer_solution.get_viscosity(temp, shear)
        else:
            return self._inferred_viscosity

    def set_inferred_viscosity(self, v: float):
        self._inferred_viscosity = v

    def get_density(self, temp: float):
        # print('filler density at ', temp, 'deg C')
        return self.oil_sample.ref_objects[0].get_density(temp)


class SpecificFluidView(QDialog, utils.SignalView):

    def __init__(self, parent, injection_fluid: InjectionFluid=None):
        super(SpecificFluidView, self).__init__(parent=parent)

        self.rheology_tool = None
        self.phri_tool = None
        self.mixing_sheet_tool = None

        mode = 'build'
        if injection_fluid is not None:
            mode = 'view'

        self.injection_fluid = injection_fluid

        # sf = parent.width() / 280.
        try:
            sf = parent.parent().sf
        except Exception as e:
            print(e)
            sf = 1.
        self.sf = sf
        self.setFixedWidth(int(sf * 350))
        self.setFixedHeight(int(sf * 130))
        layout = QVBoxLayout()
        self.setLayout(layout)

        self.name_text = QLabel(parent=self, text='Injection Fluid Name:')
        layout.addWidget(self.name_text)
        self.fluid_name_edit = QLineEdit(parent=self)
        layout.addWidget(self.fluid_name_edit)

        combo_box = QComboBox(parent=self)
        layout.addWidget(combo_box)

        self.add_button = QPushButton(parent=parent, text='+')
        self.add_button.setFixedSize(int(sf * 25), int(sf * 25))
        self.add_button.clicked.connect(self.add_additive)
        self.remove_button = QPushButton(parent=parent, text='-')
        self.remove_button.setFixedSize(int(sf * 25), int(sf * 25))
        self.remove_button.clicked.connect(self.remove_additive)
        self.submit_button = QPushButton(parent=parent, text='Submit')
        self.submit_button.setFixedHeight(int(sf * 25))
        self.submit_button.clicked.connect(self.submit)
        self.submit_button.clicked.connect(self.close)
        h_layout = QHBoxLayout()
        layout.addLayout(h_layout)
        h_layout.addWidget(self.add_button)
        h_layout.addWidget(self.remove_button)
        h_layout.addWidget(self.submit_button)

        self.source_fluid_combo_box = combo_box
        self.combo_widgets = []

        if mode == 'view':

            self.fluid_name_edit.setReadOnly(True)
            self.add_button.setEnabled(False)
            self.remove_button.setEnabled(False)
            self.submit_button.setEnabled(False)
            self.source_fluid_combo_box.setEnabled(False)

    def remove_additive(self):

        if self.combo_widgets:
            w = self.combo_widgets.pop()
            w.close()
            self.layout().removeWidget(w)
            del w
            self.setFixedHeight(self.height() - int(self.sf * 51))

    def launch_mixing_sheet_tool(self):
        if isinstance(self.injection_fluid.specific_fluid, BrineInjectionFluid):
            self.mixing_sheet_tool = MixingSheetTool(self.parent(), self.injection_fluid.specific_fluid)
            self.mixing_sheet_tool.close_signal.connect(self.mixing_sheet_tool_closed)

    def mixing_sheet_tool_closed(self, _):
        self.mixing_sheet_tool = None

    def launch_phri_tool(self):
        if isinstance(self.injection_fluid.specific_fluid, BrineInjectionFluid):
            self.phri_tool = PHRITool(self.parent(), self.injection_fluid.specific_fluid)
            self.phri_tool.close_signal.connect(self.phri_tool_closed)

    def phri_tool_closed(self, _):
        self.phri_tool = None

    def rheology_tool_closed(self, _):
        self.rheology_tool = None

    def closeEvent(self, a0) -> None:
        super(SpecificFluidView, self).closeEvent(a0)
        self.close_signal.emit(self)
        if self.rheology_tool is not None:
            self.rheology_tool.close()
        if self.phri_tool is not None:
            self.phri_tool.close()
        if self.mixing_sheet_tool is not None:
            self.mixing_sheet_tool.close()


class BrineInjectionFluidView(SpecificFluidView):

    def __init__(self, parent, brine_if: InjectionFluid = None):
        super(BrineInjectionFluidView, self).__init__(parent=parent, injection_fluid=brine_if)

        # this is probably built on an incorrect premise...the "else" likely executes all the time, and works in the
        # subsequent code
        print('BrineInjectionFluidView parent:', parent.parent().parent())
        if isinstance(self.parent(), PolymerTool):
            self.pml = self.parent().parent()
        elif isinstance(self.parent().parent(), PolymerTool):
            self.pml = self.parent().parent().parent()
        else:
            self.pml = self.parent().parent().parent().parent()

        source_names = self.pml.brine_list.signal_lists[0].item_names
        self.source_fluid_combo_box.addItems(source_names)

        mode = 'build'
        if brine_if is not None:
            mode = 'view'

        self.setWindowTitle('Brine-based Injection Fluid')

        if mode == 'view':

            brine_objects = self.pml.brine_list.signal_lists[0].objects
            for i, brine_object in enumerate(brine_objects):
                if brine_object == brine_if.ref_objects[0]:
                    self.source_fluid_combo_box.setCurrentIndex(i)
            if len(brine_if.ref_objects) > 1:
                additives = brine_if.ref_objects[1:]
                for i, additive in enumerate(additives):
                    self.add_additive(i)
                    w = self.combo_widgets[-1]
                    w.blockSignals(True)

                    if isinstance(additive, Polymer):
                        index = 0
                    elif isinstance(additive, CoSolvent):
                        index = 1
                    elif isinstance(additive, Surfactant):
                        index = 2
                    elif isinstance(additive, Scavenger):
                        index = 3
                    elif isinstance(additive, Other):
                        index = 4
                    elif isinstance(additive, Formulation):
                        index = 5

                    w.combo_box.setCurrentIndex(index)
                    additive_stocks = w.source_lists[0].signal_lists[index].objects
                    for j, additive_stock in enumerate(additive_stocks):
                        if additive_stock == additive:
                            w.combo_box_2.setCurrentIndex(j+1)
                            w.edit.setText(str(int(brine_if.concentrations[i])))

                    w.blockSignals(False)

            if isinstance(brine_if.specific_fluid, BrineInjectionFluid):
                if not isinstance(self.parent().parent(), PolymerTool):
                    try:
                        self.launch_mixing_sheet_tool()
                        self.launch_phri_tool()
                    except Exception as e:
                        QMessageBox(parent=self.parent(), text=str(e)).exec_()
                self.fluid_name_edit.setText(brine_if.name)
                ps = brine_if.specific_fluid.polymer_solution
                if ps is not None:
                    if isinstance(self.parent().parent(), PolymerTool):
                        PolymerRheologyImportTool(ps, self.parent())
                    else:
                        try:
                            ps.primary_additive.estimate_viscosity(ps.base_fluid, 75., 7.3, ps.concentration)
                            self.rheology_tool = PolymerSolutionView(brine_if.specific_fluid.polymer_solution, self.parent())
                            self.rheology_tool.close_signal.connect(self.rheology_tool_closed)
                        except Exception as e:
                            QMessageBox(parent=self, text=str(e)).exec_()

        self.show()

        if isinstance(parent.parent(), PolymerTool) and mode == 'build':
            self.add_additive()
            self.combo_widgets[0].combo_box.setEnabled(False)
            self.combo_widgets[0].combo_box_2.setCurrentIndex(self.pml.chemical_list.signal_lists[0].current_row + 1)
            self.combo_widgets[0].combo_box_2.setEnabled(False)

    def add_additive(self, conc_index: int=-1):
        new_widget = AdditiveWidget([self.pml.chemical_list], self, conc_index, True)
        self.layout().addWidget(new_widget)
        if self.injection_fluid is not None:
            new_widget.edit.editingFinished.connect(new_widget.update_injection_fluid_concentration)
            new_widget.combo_box.setEnabled(False)
            new_widget.combo_box_2.setEnabled(False)
        self.combo_widgets.append(new_widget)
        self.setFixedHeight(self.height() + int(self.sf * 51))

    def submit(self):

        name = self.fluid_name_edit.text()
        if not name:
            return

        k = self.source_fluid_combo_box.currentIndex()
        brine_obj = self.pml.brine_list.signal_lists[0].objects[k]

        concentrations = []
        ref_objs = [brine_obj]
        for w in self.combo_widgets:
            i = w.combo_box.currentIndex()
            j = w.combo_box_2.currentIndex()
            if j > 0:
                ref_obj = w.source_lists[0].signal_lists[i].objects[j-1]
                ref_objs.append(ref_obj)
            concentrations.append(float(w.edit.text()))

        l_before = len(self.parent().pm_list.signal_lists[0].objects)
        self.parent().submit_name_ref_objects(name, ref_objs, [concentrations, 'brine-based'])
        obj_after = self.parent().pm_list.signal_lists[0].objects

        if obj_after:
            # for obj in obj_after:
            #     print(obj, obj.specific_fluid.polymer_solution)
            if len(obj_after) == l_before + 1:
                obj = obj_after[-1]
                # print('parent for injection fluid view: ', self.parent().parent())
                InjectionFluidView(obj, self.parent().parent())
                # if obj.specific_fluid.polymer_solution is not None:
                #     PolymerSolutionView(obj.specific_fluid.polymer_solution, self.parent())


class OilInjectionFluidView(SpecificFluidView):

    def __init__(self, parent, oil_if: InjectionFluid = None):
        super(OilInjectionFluidView, self).__init__(parent=parent, injection_fluid=oil_if)

        if isinstance(parent, OilSampleTool):
            print('first case executing in OilInjectionFluidView')
            pmv = parent.parent()
            print(pmv)
            is_true_if = False
        elif isinstance(parent.parent(), OilSampleTool):
            print('second case executing in OilInjectionFluidView')
            pmv = parent.parent().parent()
            is_true_if = False
        else:
            print('third case executing in OilInjectionFluidView')
            pmv = parent.parent().parent().parent()
            is_true_if = True

        # if hasattr(parent.parent().parent().parent(), 'oil_list'):
        #     pmv = parent.parent().parent().parent()
        #     is_true_if = True
        # elif hasattr(parent.parent().parent(), 'oil_list'):
        #     pmv = parent.parent().parent()
        #     is_true_if = False
        # else:
        #     print('third case executing in OilInjectionFluidView')
        #     pmv = parent.parent()
        #     print(pmv)
        #     is_true_if = False

        self.pmv = pmv
        self.is_true_if = is_true_if
        self.enable_combo_box = True

        source_names = pmv.oil_list.signal_lists[1].item_names
        self.source_fluid_combo_box.addItems(source_names)

        mode = 'build'
        if oil_if is not None:
            mode = 'view'
            self.fluid_type = oil_if.get_fluid_type()
        else:
            self.add_button.setEnabled(False)
            self.remove_button.setEnabled(False)
            self.fluid_type = FluidType.OIL

        if is_true_if or oil_if is None:
            self.setWindowTitle('Oil-based Injection Fluid')
        else:
            self.setWindowTitle(oil_if.specific_fluid.oil_sample.name)
        self.slis = [2, 0]

        if mode == 'view':

            oil_objects = pmv.oil_list.signal_lists[1].objects
            for i, oil_object in enumerate(oil_objects):
                if oil_object == oil_if.ref_objects[0]:
                    self.source_fluid_combo_box.setCurrentIndex(i)
            if len(oil_if.ref_objects) > 1:
                additives = oil_if.ref_objects[1:]
                for i, additive in enumerate(additives):
                    self.add_additive(i)
                    w = self.combo_widgets[-1]
                    w.blockSignals(True)

                    index = 1

                    if isinstance(additive, Diluent):
                        index = 0
                    elif isinstance(additive, Gas):
                        index = 1

                    w.combo_box.setCurrentIndex(index)
                    additive_stocks = w.source_lists[index].signal_lists[self.slis[index]].objects
                    for j, additive_stock in enumerate(additive_stocks):
                        if additive_stock == additive:
                            w.combo_box_2.setCurrentIndex(j + 1)
                            w.edit.setText(str(int(oil_if.concentrations[i])))

                    w.blockSignals(False)

            if isinstance(oil_if.specific_fluid, OilInjectionFluid):
                self.fluid_name_edit.setText(oil_if.name)
                if oil_if.specific_fluid.polymer_solution is not None:
                    psv = PolymerSolutionView(oil_if.specific_fluid.polymer_solution, self.parent())
                    self.rheology_tool = psv
                    self.rheology_tool.close_signal.connect(self.rheology_tool_closed)
                    if is_true_if and not oil_if.specific_fluid.additives:
                        psv.submit_button.clicked.disconnect()
                        psv.allow_data_exclusion = False
                        rheo_fluids = oil_if.specific_fluid.oil_sample.rheologies_list.signal_lists[0].objects
                        if rheo_fluids:
                            psv.tool_parent = self
                            psv.oil_if = oil_if
                            psv.submit_button.clicked.connect(psv.rheo_selection_tool_launcher)

                    psv.averaging_checkbox.setChecked(True)
                    psv.toggle_averaging(True)
                    psv.averaging_checkbox.setEnabled(False)
                    # print('averaging on (to start)?:', oil_if.specific_fluid.polymer_solution.averaging_on)
                else:
                    self.setWindowTitle('Live Oil Injection Fluid')

        if mode == 'view' and isinstance(parent.parent(), OilSampleTool):
            self.close()
        else:
            self.show()

    def add_additive(self, conc_index: int=-1):
        new_widget = AdditiveWidget([self.pmv.oil_list,
                                     self.pmv.gas_list], self, conc_index, False)
        new_widget.combo_box.setEnabled(self.enable_combo_box)
        self.layout().addWidget(new_widget)
        if self.injection_fluid is not None:
            new_widget.edit.editingFinished.connect(new_widget.update_injection_fluid_concentration)
            new_widget.combo_box.setEnabled(False)
            new_widget.combo_box_2.setEnabled(False)
        self.combo_widgets.append(new_widget)
        self.setFixedHeight(self.height() + int(self.sf * 51))

    def un_pack_from_source_oil(self, oil: Oil):

        gases, ppms = oil.get_live_oil_composition()

        for gas, ppm in zip(gases, ppms):
            self.add_additive()
            widget = self.combo_widgets[-1]
            widget.combo_box.setCurrentIndex(1)
            widget.combo_box.setEnabled(False)

    def submit_wrapper(self):
        try:
            self.submit()
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()

    def submit(self):

        name = self.fluid_name_edit.text()
        if not name:
            return

        k = self.source_fluid_combo_box.currentIndex()
        if isinstance(self.parent(), OilSampleTool):
            oil_obj = self.parent().parent().oil_list.signal_lists[1].objects[k]
        elif isinstance(self.parent().parent(), OilSampleTool):
            oil_obj = self.parent().parent().parent().oil_list.signal_lists[1].objects[k]
        else:
            oil_obj = self.parent().parent().parent().parent().oil_list.signal_lists[1].objects[k]

        concentrations = []
        additives = []
        ref_objs = [oil_obj]

        if self.fluid_type == FluidType.LIVE_OIL:
            gases, ppms = oil_obj.ref_objects[0].get_live_oil_composition()
            if not gases:
                QMessageBox(parent=self, text='No gas composition for live oil.').exec_()
                return
            for gas, ppm in zip(gases, ppms):
                if isnan(ppm):
                    QMessageBox(parent=self, text='NaN gas concentration encountered.\n'
                                                  'Check oil API gravity, GOR, and gas gravity').exec_()
                    return
                concentrations.append(ppm)
                additives.append(gas)
                ref_objs.append(gas)

        for w in self.combo_widgets:
            i = w.combo_box.currentIndex()
            j = w.combo_box_2.currentIndex()
            if j > 0:
                ref_obj = w.source_lists[i].signal_lists[self.slis[i]].objects[j-1]
                ref_objs.append(ref_obj)
                additives.append(ref_obj)
            concentrations.append(float(w.edit.text()))

        if not self.is_true_if and not oil_obj.composition_set:
            oil_obj.additives = additives
            oil_obj.concentrations = concentrations
            oil_obj.composition_set = True
            print('viewing composition from OilSampleTool')
            self.parent().view_composition()
            print('completed')
            return

        l_before = len(self.parent().pm_list.signal_lists[0].objects)
        # print('OilInjectionFluid submit: ', name, ref_objs, [concentrations, 'oil-based'])
        self.parent().submit_name_ref_objects(name, ref_objs, [concentrations, 'oil-based'])
        obj_after = self.parent().pm_list.signal_lists[0].objects

        if obj_after:
            # for obj in obj_after:
            #     print(obj, obj.specific_fluid.polymer_solution)
            if len(obj_after) == l_before + 1:
                obj = obj_after[-1]
                # print('parent for injection fluid view: ', self.parent().parent())
                try:
                    InjectionFluidView(obj, self.parent().parent())
                except Exception as e:
                    QMessageBox(parent=self, text=str(e)).exec_()
                # if obj.specific_fluid.polymer_solution is not None:
                #     PolymerSolutionView(obj.specific_fluid.polymer_solution, self.parent())


class AdditiveWidget(QWidget):

    def __init__(self, source_lists: list, parent=None, conc_index=-1, is_brine: bool=True):
        super(AdditiveWidget, self).__init__(parent=parent)

        self.source_lists = source_lists
        self.conc_index = conc_index
        self.is_brine = is_brine

        h_layout = QHBoxLayout()
        h_layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(h_layout)
        self.combo_box = QComboBox(parent=self)
        if is_brine:
            self.combo_box.addItems(['Polymer', 'Co-solvent', 'Surfactant', 'Scavenger', 'Other', 'Formulation'])
        else:
            self.combo_box.addItems(['Diluent', 'Gas'])
            self.slis = [2, 0]
        self.combo_box.currentIndexChanged.connect(self.set_combobox_names)
        self.combo_box_2 = QComboBox(parent=self)
        if self.is_brine:
            self.combo_box_2.currentIndexChanged.connect(self.set_formulation_concentration)

        self.edit = QLineEdit(parent=self)
        self.edit.setFixedWidth(int(parent.sf * 60))
        self.edit.editingFinished.connect(self.check_edit_value)
        text = QLabel(parent=self, text='[ppm]')
        text.setFixedWidth(int(parent.sf * 40))
        h_layout.addWidget(self.combo_box)
        h_layout.addWidget(self.combo_box_2)
        h_layout.addWidget(self.edit)
        h_layout.addWidget(text)
        self.concentration_text = text

        self.set_combobox_names()

    def check_edit_value(self):

        txt = self.edit.text()

        try:
            val = float(txt)
            if val < 0.:
                raise ValueError
        except ValueError:
            self.edit.setText('')

    def set_combobox_names(self):

        self.combo_box_2.blockSignals(True)
        self.combo_box_2.clear()
        cb_list = list(' ')
        i = self.combo_box.currentIndex()
        if self.is_brine:
            names = self.source_lists[0].signal_lists[i].item_names
        else:
            names = self.source_lists[i].signal_lists[self.slis[i]].item_names
        for name in names:
            cb_list.append(name)

        self.combo_box_2.addItems(cb_list)
        self.combo_box_2.blockSignals(False)

    def set_combobox_to_name(self, name: str):

        i = self.combo_box_2.findText(name)
        if i != -1:
            self.combo_box_2.setCurrentIndex(i)

    def set_formulation_concentration(self):

        if self.combo_box.currentIndex() != 5:
            self.edit.setEnabled(True)
            return

        j = self.combo_box_2.currentIndex() - 1
        formulation = self.source_lists[0].signal_lists[5].objects[j]
        val = np_sum(formulation.concentrations)
        self.edit.setText(str(val))
        self.edit.setEnabled(False)

    def update_injection_fluid_concentration(self):

        if self.parent().injection_fluid is None:
            return

        if self.conc_index == -1:
            return

        # print(self.conc_index, float(self.edit.text()))
        # print('concentrations before: ', self.parent().injection_fluid.concentrations)
        self.parent().injection_fluid.concentrations[self.conc_index] = float(self.edit.text())
        if self.combo_box.currentIndex() == 0 and self.is_brine:
            self.parent().injection_fluid.specific_fluid.polymer_solution.concentration = float(self.edit.text())
        # print('concentrations after: ', self.parent().injection_fluid.concentrations)


class AdditiveWidgetGasMix(AdditiveWidget):

    def __init__(self, source_lists: list, parent=None, conc_index=-1):
        super(AdditiveWidgetGasMix, self).__init__(source_lists=source_lists, parent=parent, conc_index=conc_index,
                                                   is_brine=False)
        self.combo_box.setCurrentIndex(1)
        self.combo_box.setEnabled(False)
        self.concentration_text.setText('[%wt]')


class PHRITool(QDialog, utils.SignalView):

    def __init__(self, parent, brine_if: BrineInjectionFluid):
        super(PHRITool, self).__init__(parent=parent)
        sf = parent.parent().sf
        self.sf = sf
        self.brine_if = brine_if

        self.setWindowTitle(brine_if.injection_fluid.name + ' pH/R.I. Tool')
        self.setFixedSize(int(sf * 250), int(sf * 75))
        lyt = QGridLayout()
        self.setLayout(lyt)

        ph_label = QLabel(parent=self, text='pH:')
        ri_label = QLabel(parent=self, text='R.I. [o/oo]:')
        self.ph_edit = QLineEdit(parent=self)
        self.ph_edit.editingFinished.connect(lambda: self.vet_edit(self.ph_edit))
        self.ri_edit = QLineEdit(parent=self)
        self.ri_edit.editingFinished.connect(lambda: self.vet_edit(self.ri_edit))

        lyt.addWidget(ph_label, 0, 0, 1, 1)
        lyt.addWidget(ri_label, 1, 0, 1, 1)
        lyt.addWidget(self.ph_edit, 0, 1, 1, 1)
        lyt.addWidget(self.ri_edit, 1, 1, 1, 1)

        self.sync_edits_to_if()
        self.show()

    def vet_edit(self, edit: QLineEdit):

        txt = edit.text()

        try:
            val = float(txt)
            if val < 0.:
                raise ValueError

            if edit == self.ph_edit:
                self.brine_if.pH = val
            elif edit == self.ri_edit:
                self.brine_if.RI = val

        except ValueError:
            self.sync_edits_to_if()

    def sync_edits_to_if(self):

        self.ph_edit.clear()
        self.ri_edit.clear()

        if not isnan(self.brine_if.pH):
            self.ph_edit.setText(str(self.brine_if.pH))

        if not isnan(self.brine_if.RI):
            self.ri_edit.setText(str(self.brine_if.RI))

    def closeEvent(self, a0) -> None:
        super(PHRITool, self).closeEvent(a0)
        self.close_signal.emit(self)


class MixingSheetTool(QDialog, utils.SignalView):

    def __init__(self, parent, brine_if: BrineInjectionFluid):
        super(MixingSheetTool, self).__init__(parent=parent)
        sf = parent.parent().sf
        self.sf = sf
        self.brine_if = brine_if
        self.stocks_mixing_list = {}
        self.fluids_mixing_list = []
        self.fluids_mixing_salinity = 0.
        self.balance_brine = []
        self.setLayout(QVBoxLayout())
        self.setWindowTitle(brine_if.injection_fluid.name + ' Mixing Sheet Tool')

        self.resize(int(415 * sf), int(625 * sf))

        pmv = parent.parent().parent().parent()
        self.available_brines = pmv.brine_list.signal_lists[0].objects
        self.available_brines_names = []
        for brine in self.available_brines:
            self.available_brines_names.append(brine.name)

        self.mixing_stocks_panel = None
        if brine_if.has_formulation() or brine_if.has_additive_salt() or brine_if.has_polymer():
            self.mixing_stocks_panel = MixingStocksPanel(parent=self, brine_if=brine_if)
            stocks_label = QLabel(parent=self, text='Stocks')
            label_font = stocks_label.font()
            label_font.setPointSize(10)
            label_font.setBold(True)
            stocks_label.setFont(label_font)
            stocks_label.setFixedHeight(18)
            self.layout().addWidget(stocks_label)
            self.layout().addWidget(self.mixing_stocks_panel)

        fluid_mixing_label = QLabel(parent=self, text=brine_if.injection_fluid.name + ' Fluid Mixing')
        label_font = fluid_mixing_label.font()
        label_font.setPointSize(10)
        label_font.setBold(True)
        fluid_mixing_label.setFont(label_font)
        fluid_mixing_label.setFixedHeight(18)
        self.layout().addWidget(fluid_mixing_label)
        self.fluid_mixing_panel = MixingFluidsPanel(parent=self, brine_if=brine_if)
        self.layout().addWidget(self.fluid_mixing_panel)

        if self.mixing_stocks_panel is not None:
            self.mixing_stocks_panel.stockUpdated.connect(self.fluid_mixing_panel.update_final_mixing)
            print(self.mixing_stocks_panel, self.mixing_stocks_panel.get_my_natural_height(), self.mixing_stocks_panel.height())

        pdf_button = QPushButton(parent=self, text='Generate Mixing Sheet')
        pdf_button.clicked.connect(self.generate_mixing_sheet)
        self.layout().addWidget(pdf_button)
        self.mixing_sheet_button = pdf_button

        self.show()
        self.resize(self.width(), self.get_my_natural_height())
        self.setFixedWidth(self.width())
        print('Mixing Sheet Tool height is {}.'.format(self.get_my_natural_height()))
        for stock, stock_c in brine_if.stocks_list.items():
            print('{}: {}'.format(stock, stock_c))

    def get_my_natural_height(self) -> int:

        height = self.fluid_mixing_panel.get_my_natural_height()
        print('fluid mixing panel height = {}'.format(height))
        if self.mixing_stocks_panel is not None:
            print('mixing stocks panel height = {}'.format(self.mixing_stocks_panel.get_my_natural_height()))
            height += self.mixing_stocks_panel.get_my_natural_height()
        return int(height) + 30

    def generate_mixing_sheet(self):
        file_name = utils.file_open_dlg(self.brine_if.injection_fluid.name + ' Mixing Sheet', '(PDF) *.pdf', True)

        try:
            snl = []
            sml = []
            for key, value in self.stocks_mixing_list.items():
                snl.append(key)
                sml.append(value)
            reports.create_injection_fluid_mixing_sheet(snl, sml, self.balance_brine,
                                                        self.brine_if.injection_fluid.name, self.fluids_mixing_list,
                                                        self.fluids_mixing_salinity, file_name)
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()

    def closeEvent(self, a0) -> None:
        super(MixingSheetTool, self).closeEvent(a0)
        self.close_signal.emit(self)


class MixingValueEdit(QLineEdit):

    valueChanged = pyqtSignal(object)

    def __init__(self, parent, key, val: float=0., check_val: float=0.):
        super(MixingValueEdit, self).__init__(parent=parent)

        self.check_val = check_val
        self.val = val
        self.setText(str(val))
        self.check_value()
        self.editingFinished.connect(self.check_value)
        self.key = key

    def check_value(self):

        txt = self.text()
        try:
            new_val = float(txt)
            if new_val < self.check_val:
                raise ValueError
            if new_val != self.val:
                self.val = new_val
                self.valueChanged.emit(self)

        except ValueError:
            self.setText(str(self.check_val))


class MixingCombo(QComboBox):

    valueChanged = pyqtSignal(object)

    def __init__(self, parent, key=None):
        super(MixingCombo, self).__init__(parent=parent)
        self.currentTextChanged.connect(self.text_changed)
        self.key = key

    def text_changed(self, _):
        self.valueChanged.emit(self)


class MixingStocksPanel(QWidget):

    stockUpdated = pyqtSignal()

    def __init__(self, parent, brine_if: BrineInjectionFluid):
        super(MixingStocksPanel, self).__init__(parent=parent)
        lyt = QGridLayout()
        self.setLayout(QVBoxLayout())
        self.frame = QFrame(parent=self)
        self.frame.setFrameShape(QFrame.StyledPanel)
        self.frame.setLineWidth(1)
        self.frame.setLayout(lyt)
        self.n_rows = 0.
        self.row_height = 0.
        self.layout().addWidget(self.frame)

        self.formulation_concentration_labels = []
        self.formulation_concentration_edits = []
        self.formulation_base_labels = []
        self.formulation_base_combos = []
        self.formulation_component_concentration_edits = []
        self.formulation_mass_edits = []
        self.formulation_mixing_labels = []
        self.salt_concentration_labels = []
        self.salt_concentration_edits = []
        self.salt_base_labels = []
        self.salt_base_combos = []
        self.polymer_concentration_labels = []
        self.polymer_concentration_edits = []
        self.polymer_base_labels = []
        self.polymer_base_combos = []

        i, nearest_brine = brine_if.brine.find_nearest_base_composition_brine(parent.available_brines)
        row = 0

        if brine_if.has_formulation():
            j = 0
            formulations, concentrations = brine_if.formulations_and_concentrations()
            for formulation, c in zip(formulations, concentrations):
                c_label = QLabel(parent=self, text=formulation.name + ' Conc. (X):')
                val = brine_if.stocks_list[formulation]
                c_edit = MixingValueEdit(parent=self, key=formulation, val=val, check_val=1.)
                self.row_height = c_edit.height()
                c_edit.valueChanged.connect(self.update_formulation_mixing_wrapper)
                b_label = QLabel(parent=self, text='Base Brine:')
                b_combo = MixingCombo(parent=self, key=formulation)
                b_combo.addItems(parent.available_brines_names)
                try:
                    k = parent.available_brines.index(brine_if.stocks_base_fluids[formulation])
                except Exception as e:
                    print(e)
                    k = i
                    brine_if.stocks_base_fluids[formulation] = nearest_brine

                b_combo.setCurrentIndex(k)
                b_combo.valueChanged.connect(self.update_formulation_mixing_wrapper_combo)
                m_label = QLabel(parent=self, text='Total Mass [g]:')
                val = brine_if.formulation_mixing_masses[formulation]
                m_edit = MixingValueEdit(parent=self, key=formulation, val=val, check_val=0.)
                m_edit.valueChanged.connect(self.update_formulation_mixing_wrapper)

                self.formulation_concentration_labels.append(c_label)
                self.formulation_concentration_edits.append(c_edit)
                self.formulation_base_labels.append(b_label)
                self.formulation_base_combos.append(b_combo)
                self.formulation_mass_edits.append(m_edit)

                lyt.addWidget(c_label, row, 0, 1, 1)
                lyt.addWidget(c_edit, row, 1, 1, 1)
                row += 1
                lyt.addWidget(b_label, row, 0, 1, 1)
                lyt.addWidget(b_combo, row, 1, 1, 1)
                row += 1

                edit_list = []
                for item in formulation.sc_list:
                    cc_label = QLabel(parent=self, text=item.name + ' Stock Conc. [%]')
                    val = brine_if.stocks_list[item]
                    cc_edit = MixingValueEdit(parent=self, key=item, val=val, check_val=0.)
                    cc_edit.valueChanged.connect(self.update_formulation_mixing_wrapper)
                    edit_list.append(cc_edit)
                    lyt.addWidget(cc_label, row, 0, 1, 1)
                    lyt.addWidget(cc_edit, row, 1, 1, 1)
                    row += 1

                self.formulation_component_concentration_edits.append(edit_list)

                lyt.addWidget(m_label, row, 0, 1, 1)
                lyt.addWidget(m_edit, row, 1, 1, 1)
                row += 1

                mix_label = QLabel(parent=self, text='')
                font = mix_label.font()
                font.setFamily('Courier')
                # font.setPointSize(12)
                mix_label.setFont(font)
                self.formulation_mixing_labels.append(mix_label)
                lyt.addWidget(mix_label, row, 0, 1, 2)
                row += 1
                self.update_formulation_mixing(j)
                j += 1

        if brine_if.has_additive_salt():
            salts, concentrations = brine_if.brine.additive_salts_and_concentrations()

            for salt, c in zip(salts, concentrations):

                if salt not in brine_if.stocks_list.keys():
                    if c < 10.:
                        brine_if.stocks_list[salt] = 10.
                    else:
                        brine_if.stocks_list[salt] = c

                c_label = QLabel(parent=self, text=salt + ' Conc. [wt%]:')
                val = brine_if.stocks_list[salt]
                c_edit = MixingValueEdit(parent=self, key=salt, val=val, check_val=c)
                self.row_height = c_edit.height()
                c_edit.valueChanged.connect(self.generate_update_signal)
                b_label = QLabel(parent=self, text='Base Brine:')
                b_combo = MixingCombo(parent=self, key=salt)
                b_combo.addItems(parent.available_brines_names)
                try:
                    k = parent.available_brines.index(brine_if.stocks_base_fluids[salt])
                except Exception as e:
                    print(e)
                    k = i
                    brine_if.stocks_base_fluids[salt] = nearest_brine

                b_combo.setCurrentIndex(k)
                b_combo.valueChanged.connect(self.generate_update_signal)
                self.salt_concentration_labels.append(c_label)
                self.salt_concentration_edits.append(c_edit)
                self.salt_base_labels.append(b_label)
                self.salt_base_combos.append(b_combo)
                lyt.addWidget(c_label, row, 0, 1, 1)
                lyt.addWidget(c_edit, row, 1, 1, 1)
                row += 1
                lyt.addWidget(b_label, row, 0, 1, 1)
                lyt.addWidget(b_combo, row, 1, 1, 1)
                row += 1

        if brine_if.has_polymer():
            polymers, concentrations = brine_if.polymers_and_concentrations()
            for polymer, c in zip(polymers, concentrations):
                c_label = QLabel(parent=self, text=polymer.name + ' Conc. [ppm]:')
                val = brine_if.stocks_list[polymer]
                c_edit = MixingValueEdit(parent=self, key=polymer, val=val, check_val=c)
                self.row_height = c_edit.height()
                c_edit.valueChanged.connect(self.generate_update_signal)
                b_label = QLabel(parent=self, text='Base Brine:')
                b_combo = MixingCombo(parent=self, key=polymer)
                b_combo.addItems(parent.available_brines_names)
                try:
                    k = parent.available_brines.index(brine_if.stocks_base_fluids[polymer])
                except Exception as e:
                    print(e)
                    k = i
                    brine_if.stocks_base_fluids[polymer] = nearest_brine

                b_combo.setCurrentIndex(k)
                b_combo.valueChanged.connect(self.generate_update_signal)
                self.polymer_concentration_labels.append(c_label)
                self.polymer_concentration_edits.append(c_edit)
                self.polymer_base_labels.append(b_label)
                self.polymer_base_combos.append(b_combo)
                lyt.addWidget(c_label, row, 0, 1, 1)
                lyt.addWidget(c_edit, row, 1, 1, 1)
                row += 1
                lyt.addWidget(b_label, row, 0, 1, 1)
                lyt.addWidget(b_combo, row, 1, 1, 1)
                row += 1

        self.n_rows = row + 1

    def get_my_natural_height(self):

        return (self.n_rows + len(self.formulation_mixing_labels)) * self.row_height

    def check_formulation_concentration(self):

        txt = self.formulation_concentration_edit.text()
        try:
            val = float(txt)
            if val < 1.:
                raise ValueError

        except ValueError:
            self.formulation_concentration_edit.setText('8')

    def find_edit(self, edit: MixingValueEdit) -> int:

        if edit in self.formulation_concentration_edits:
            i = self.formulation_concentration_edits.index(edit)
        elif edit in self.formulation_mass_edits:
            i = self.formulation_mass_edits.index(edit)
        else:
            i = 0
            for j, group in enumerate(self.formulation_component_concentration_edits):
                if edit in group:
                    i = j

        print('c_edit index is {}'.format(i))
        return i

    def find_combo(self, combo: MixingCombo) -> int:

        return self.formulation_base_combos.index(combo)

    def check_edit_value(self, c_edit: MixingValueEdit, c: float):

        txt = c_edit.text()
        try:
            val = float(txt)
            if val < c:
                raise ValueError

            if c_edit in [*self.formulation_component_concentration_edits, *self.formulation_concentration_edits,
                          *self.formulation_mass_edits]:

                if c_edit in self.formulation_concentration_edits:
                    i = self.formulation_concentration_edits.index(c_edit)
                elif c_edit in self.formulation_mass_edits:
                    i = self.formulation_mass_edits.index(c_edit)
                else:
                    i = 0
                    for j, group in enumerate(self.formulation_component_concentration_edits):
                        if c_edit in group:
                            i = j

                print('c_edit index is {}'.format(i))

                self.update_formulation_mixing(i)

        except ValueError:
            c_edit.setText(str(c))

    def update_formulation_mixing_wrapper(self, edit: MixingValueEdit):

        self.update_formulation_mixing(self.find_edit(edit=edit))
        self.generate_update_signal(obj=edit)

    def update_formulation_mixing_wrapper_combo(self, combo: MixingCombo):

        self.update_formulation_mixing(self.find_combo(combo=combo))
        self.generate_update_signal(obj=combo)

    def update_formulation_mixing(self, i: int):

        label_text = ''
        label = self.formulation_mixing_labels[i]
        stocks_mixing_list = []
        cc_edits = self.formulation_component_concentration_edits[i]
        formulations, _ = self.parent().brine_if.formulations_and_concentrations()
        formulation = formulations[i]
        concentration_factor = float(self.formulation_concentration_edits[i].text())
        total_mass = float(self.formulation_mass_edits[i].text())
        order = 0
        while total_mass >= 10. ** (order + 1):
            order += 1
        brine_mass = total_mass

        max_len = 0.
        for item, cc_edit in zip(formulation.sc_list, cc_edits):
            name_text = '{} ({:.1f}%)'.format(item.name, float(cc_edit.text()))
            if len(name_text) > max_len:
                max_len = len(name_text)
        if max_len < len(self.formulation_base_combos[i].currentText()):
            max_len = len(self.formulation_base_combos[i].currentText())

        for item, c, cc_edit in zip(formulation.sc_list, formulation.concentrations, cc_edits):
            c_factor = 100. / float(cc_edit.text())
            stock_c = 1.e-6 * c * concentration_factor * c_factor
            brine_mass -= stock_c * total_mass
            name_text = '{} ({:.1f}%)'.format(item.name, float(cc_edit.text()))
            if len(name_text) < max_len:
                for j in range(max_len - len(name_text)):
                    name_text += ' '
            val_text = '{:.3f} g'.format(stock_c * total_mass)
            true_stock_c = float(cc_edit.text()) * 1.e4
            stocks_mixing_list.append([item.name, '', true_stock_c, c * concentration_factor, stock_c * total_mass])
            val_order = 0
            while stock_c * total_mass >= 10. ** (val_order + 1):
                val_order += 1
            if val_order < order:
                for j in range(order - val_order):
                    val_text = ' ' + val_text
            label_text += name_text + ' | ' + val_text + '\n'

        name_text = self.formulation_base_combos[i].currentText()
        if len(name_text) < max_len:
            for j in range(max_len - len(name_text)):
                name_text += ' '
        val_text = '{:.3f} g'.format(brine_mass)
        stocks_mixing_list.append([name_text, '', 1000000., 1000000., brine_mass])
        val_order = 0
        while brine_mass >= 10. ** (val_order + 1):
            val_order += 1
        if val_order < order:
            for j in range(order - val_order):
                val_text = ' ' + val_text
        label_text += name_text + ' | ' + val_text
        label.setText(label_text)

        key = formulation.name + ' (' + str(concentration_factor) + 'X)'
        self.parent().stocks_mixing_list[key] = stocks_mixing_list

        # self.generate_update_signal()

    def generate_update_signal(self, obj):

        if isinstance(obj, MixingCombo) and obj.key is not None:
            print('obj, key:', obj, obj.key)
            self.parent().brine_if.stocks_base_fluids[obj.key] = self.parent().available_brines[obj.currentIndex()]

        elif isinstance(obj.key, Surfactant) or isinstance(obj.key, CoSolvent) or isinstance(obj.key, str) or \
                isinstance(obj.key, Polymer) or isinstance(obj.key, Formulation):
            if obj not in self.formulation_mass_edits:
                self.parent().brine_if.stocks_list[obj.key] = obj.val
            else:
                self.parent().brine_if.formulation_mixing_masses[obj.key] = obj.val
            print(self.parent().brine_if.stocks_list)

        self.stockUpdated.emit()


class MixingFluidsPanel(QWidget):

    def __init__(self, parent, brine_if: BrineInjectionFluid):
        super(MixingFluidsPanel, self).__init__(parent=parent)

        self.setLayout(QVBoxLayout())
        lyt = QGridLayout()
        self.frame = QFrame(parent=self)
        self.frame.setFrameShape(QFrame.StyledPanel)
        self.frame.setLineWidth(1)
        self.frame.setLayout(lyt)
        self.n_rows = 0.
        self.row_height = 0.
        self.formulation_mixing_lines = 0.
        self.layout().addWidget(self.frame)

        row = 0

        b_combo_label = QLabel(parent=self, text='Base Brine (target):')
        b_combo = MixingCombo(parent=self)
        b_combo.addItems(parent.available_brines_names)
        try:
            i = parent.available_brines.index(brine_if.mixing_base_fluid)
        except ValueError:
            i, b = brine_if.brine.find_nearest_base_composition_brine(parent.available_brines)
            brine_if.mixing_base_fluid = b

        b_combo.setCurrentIndex(i)
        b_combo.valueChanged.connect(self.update_final_mixing_wrapper)
        lyt.addWidget(b_combo_label, row, 0, 1, 1)
        lyt.addWidget(b_combo, row, 1, 1, 1)
        row += 1

        d_combo_label = QLabel(parent=self, text='Dilution Brine:')
        d_combo = MixingCombo(parent=self)
        d_combo.addItems(parent.available_brines_names)
        try:
            j = parent.available_brines.index(brine_if.mixing_dilution_fluid)
        except ValueError:
            j, b2 = brine_if.brine.find_nearest_base_composition_brine(parent.available_brines)
            brine_if.mixing_dilution_fluid = b2

        d_combo.setCurrentIndex(j)
        d_combo.valueChanged.connect(self.update_final_mixing_wrapper)
        lyt.addWidget(d_combo_label, row, 0, 1, 1)
        lyt.addWidget(d_combo, row, 1, 1, 1)
        row += 1

        d_conc_label = QLabel(parent=self, text='Dil. Brine Conc. (X)')
        val = brine_if.mixing_dilution_x_factor
        d_conc_edit = MixingValueEdit(parent=self, key=None, val=val, check_val=0.)
        self.row_height = d_conc_edit.height()
        d_conc_edit.valueChanged.connect(self.update_final_mixing_wrapper)
        lyt.addWidget(d_conc_label, row, 0, 1, 1)
        lyt.addWidget(d_conc_edit, row, 1, 1, 1)
        row += 1

        self.surfactant_stock_concentration_edits = []
        self.cosolvent_stock_concentration_edits = []

        surfactants, _ = brine_if.surfactants_and_concentrations()
        for surfactant in surfactants:
            s_label = QLabel(parent=self, text='{} Stock Conc. [wt%]: '.format(surfactant.name))
            val = brine_if.stocks_list[surfactant]
            s_edit = MixingValueEdit(parent=self, key=surfactant, val=val, check_val=0.)
            s_edit.valueChanged.connect(self.update_final_mixing_wrapper)
            lyt.addWidget(s_label, row, 0, 1, 1)
            lyt.addWidget(s_edit, row, 1, 1, 1)
            row += 1
            self.surfactant_stock_concentration_edits.append(s_edit)

        cosolvents, _ = brine_if.cosolvents_and_concentrations()
        for cosolvent in cosolvents:
            c_label = QLabel(parent=self, text='{} Stock Conc. [wt%]: '.format(cosolvent.name))
            val = brine_if.stocks_list[cosolvent]
            c_edit = MixingValueEdit(parent=self, key=cosolvent, val=val, check_val=0.)
            c_edit.valueChanged.connect(self.update_final_mixing_wrapper)
            lyt.addWidget(c_label, row, 0, 1, 1)
            lyt.addWidget(c_edit, row, 1, 1, 1)
            row += 1
            self.cosolvent_stock_concentration_edits.append(c_edit)

        mass_label = QLabel(parent=self, text='Total Mass [g]:')
        val = brine_if.mixing_mass
        mass_edit = MixingValueEdit(parent=self, key=None, val=val, check_val=0.)
        mass_edit.valueChanged.connect(self.update_final_mixing_wrapper)
        lyt.addWidget(mass_label, row, 0, 1, 1)
        lyt.addWidget(mass_edit, row, 1, 1, 1)
        row += 1

        self.n_rows = row

        mix_label = QLabel(parent=self, text='')
        font = mix_label.font()
        font.setFamily('Courier')
        mix_label.setFont(font)
        lyt.addWidget(mix_label, row, 0, 1, 2)

        self.mixing_base_combo = b_combo
        self.mixing_dilution_combo = d_combo
        self.mixing_dilution_concentration_edit = d_conc_edit
        self.mixing_mass_edit = mass_edit
        self.final_mixing_label = mix_label

        self.update_final_mixing()

    def get_my_natural_height(self):

        return (self.n_rows + self.formulation_mixing_lines) * self.row_height

    def update_final_mixing_wrapper(self, obj):

        self.update_final_mixing()

        if isinstance(obj, MixingCombo) and obj != self.mixing_base_combo and obj != self.mixing_dilution_combo:
            return

        brine_if = self.parent().brine_if

        if obj == self.mixing_base_combo:
            brine_if.mixing_base_fluid = self.parent().available_brines[obj.currentIndex()]
            return
        elif obj == self.mixing_dilution_combo:
            brine_if.mixing_dilution_fluid = self.parent().available_brines[obj.currentIndex()]
            return

        if obj == self.mixing_mass_edit:
            brine_if.mixing_mass = obj.val

        elif obj == self.mixing_dilution_concentration_edit:
            brine_if.mixing_dilution_x_factor = obj.val

        elif obj.key is not None:
            brine_if.stocks_list[obj.key] = obj.val

    def update_final_mixing(self):

        mixing_text = ''
        total_mass = self.mixing_mass_edit.val
        total_mass_order = 0
        format_string = '{:5.3f}'
        while total_mass >= 10. ** (total_mass_order + 1):
            total_mass_order += 1
            format_string = '{:' + str(int(total_mass_order + 5)) + '.3f}'
        added_mass = 0.
        brine_if = self.parent().brine_if
        brine_if.mixing_mass = self.mixing_mass_edit.val
        msp = self.parent().mixing_stocks_panel
        added_salt = 0.
        added_ions = zeros((14, ))
        brines = self.parent().available_brines
        formulation_mixing_lines = 0.
        fluids_mixing_list = []
        total_fluid_salinity = 0.

        if msp is not None:

            formulations, concentrations = brine_if.formulations_and_concentrations()
            for f, c, x_edit, b_combo in zip(formulations, concentrations, msp.formulation_concentration_edits,
                                             msp.formulation_base_combos):
                i = b_combo.currentIndex()
                brine = brines[i]
                x_factor = float(x_edit.text())
                mass = total_mass / x_factor
                added_mass += mass
                added_salt += mass * brine.get_tds()
                added_ions += mass * brine.get_composition_array()
                line = '{} | ' + format_string + ' g \n'
                mixing_text += line.format(f.name, mass)
                fluids_mixing_list.append([f.name + ' (' + brine.name + ')', '', x_factor * c, c, mass])
                formulation_mixing_lines += 1

            j = self.mixing_dilution_combo.currentIndex()
            d_salt, _ = brines[j].additive_salts_and_concentrations()
            salts, concentrations = brine_if.brine.additive_salts_and_concentrations()
            for s, c, c_edit, b_combo in zip(salts, concentrations, msp.salt_concentration_edits,
                                             msp.salt_base_combos):
                if s in d_salt:
                    continue
                i = b_combo.currentIndex()
                brine = brines[i]
                original_c = brine.add_salt_composition[s]
                stock_c = float(c_edit.text())
                brine.add_salt_composition[s] = stock_c + original_c
                mass = (c / stock_c) * total_mass
                added_mass += mass
                added_salt += mass * brine.get_tds()
                total_fluid_salinity += 10000. * c
                added_ions += mass * brine.get_composition_array()
                brine.add_salt_composition[s] = original_c
                line = '{} ({:.1f}%) | ' + format_string + ' g\n'
                mixing_text += line.format(s, stock_c, mass)
                if stock_c < 100.:
                    fluids_mixing_list.append([s + ' (' + brine.name + ')', '', stock_c * 10000., c * 10000., mass])
                else:
                    fluids_mixing_list.append([s, '', stock_c * 10000., c * 10000., mass])
                formulation_mixing_lines += 1

            polymers, concentrations = brine_if.polymers_and_concentrations()
            for p, c, c_edit, b_combo in zip(polymers, concentrations, msp.polymer_concentration_edits,
                                             msp.polymer_base_combos):
                i = b_combo.currentIndex()
                brine = brines[i]
                stock_c = float(c_edit.text())
                mass = (c / stock_c) * total_mass
                added_mass += mass
                added_salt += mass * brine.get_tds()
                added_ions += mass * brine.get_composition_array()
                line = '{} ({:.0f} ppm) | ' + format_string + ' g \n'
                mixing_text += line.format(p.name, stock_c, mass)
                fluids_mixing_list.append([p.name + ' (' + brine.name + ')', '          ', stock_c, c, mass])
                formulation_mixing_lines += 1

        scso_mass = 0.

        surfactants, concentrations = brine_if.surfactants_and_concentrations()
        for s, c, c_edit in zip(surfactants, concentrations, self.surfactant_stock_concentration_edits):
            stock_c = float(c_edit.text())
            mass = 1.e-6 * c * total_mass / (0.01 * stock_c)
            added_mass += mass
            scso_mass += mass
            line = '{} ({:.1f}%) | ' + format_string + ' g\n'
            mixing_text += line.format(s.name, stock_c, mass)
            fluids_mixing_list.append([s.name, '          ', stock_c, c, mass])
            formulation_mixing_lines += 1

        cosolvents, concentrations = brine_if.cosolvents_and_concentrations()
        for cs, c, c_edit in zip(cosolvents, concentrations, self.cosolvent_stock_concentration_edits):
            stock_c = float(c_edit.text())
            mass = 1.e-6 * c * total_mass / (0.01 * stock_c)
            added_mass += mass
            scso_mass += mass
            line = '{} ({:.1f}%) | ' + format_string + ' g\n'
            mixing_text += line.format(cs.name, stock_c, mass)
            fluids_mixing_list.append([cs.name, '          ', stock_c, c, mass])
            formulation_mixing_lines += 1

        scavengers, concentrations = brine_if.scavengers_and_concentrations()
        for s, c in zip(scavengers, concentrations):
            mass = 1.e-6 * c * total_mass
            added_mass += mass
            scso_mass += mass
            line = '{} | ' + format_string + ' g\n'
            mixing_text += line.format(s.name, mass)
            fluids_mixing_list.append([s.name, '', 1000000., c, mass])
            formulation_mixing_lines += 1

        others, concentrations = brine_if.others_and_concentrations()
        for o, c in zip(others, concentrations):
            mass = 1.e-6 * c * total_mass
            added_mass += mass
            scso_mass += mass
            line = '{} | ' + format_string + ' g\n'
            mixing_text += line.format(o.name, mass)
            fluids_mixing_list.append([o.name, '', 1000000., c, mass])
            formulation_mixing_lines += 1

        non_salts, concentrations = brine_if.brine.additive_non_salts_and_concentrations()
        for ns, c in zip(non_salts, concentrations):
            mass = 1.e-2 * c * total_mass
            added_mass += mass
            scso_mass += mass
            line = '{} | ' + format_string + ' g\n'
            mixing_text += line.format(ns, mass)
            fluids_mixing_list.append([ns, '', 1000000., 10000. * c, mass])
            formulation_mixing_lines += 1

        i = self.mixing_base_combo.currentIndex()
        brine = brines[i]
        total_fluid_salinity += brine.get_tds()
        k = self.mixing_dilution_combo.currentIndex()
        d_brine = brines[k]
        x_factor = float(self.mixing_dilution_concentration_edit.text())

        total_ions_target = brine.get_composition_array() * (total_mass - scso_mass)
        mixing_brine_ions = d_brine.get_composition_array() * x_factor
        non_zeros = where(mixing_brine_ions != 0.)
        if not list(non_zeros[0]):
            mixing_brine_mass = total_mass - added_mass
            line = '{} ({}X) | ' + format_string + ' g\n'
            mixing_text += line.format(d_brine.name, x_factor, mixing_brine_mass)
            line = 'DI | ' + format_string + ' g'
            mixing_text += line.format(0.)
            fluids_mixing_list.append([d_brine.name, '', 1000000., 1000000., mixing_brine_mass])
            fluids_mixing_list.append(['DI', '', 1000000., 1000000., 0.])
        else:
            print(total_ions_target)
            print(added_ions)
            print(mixing_brine_ions)
            masses_by_ion = divide(total_ions_target[non_zeros] - added_ions[non_zeros], mixing_brine_ions[non_zeros])
            mixing_brine_mass = np_min(masses_by_ion)
            remaining_ion_masses = total_ions_target - added_ions - mixing_brine_ions * mixing_brine_mass
            remaining_ion_composition = remaining_ion_masses / (total_mass - mixing_brine_mass - added_mass)

            print('remaining ion masses:', remaining_ion_masses)
            print('remaining ion composition:', remaining_ion_composition)
            # total_salt_target = brine.get_tds() * (total_mass - scso_mass)
            # mixing_brine_mass = (total_salt_target - added_salt) / (d_brine.get_tds() * x_factor)
            line = '{} ({}X) | ' + format_string + ' g\n'
            mixing_text += line.format(d_brine.name, x_factor, mixing_brine_mass)
            fluids_mixing_list.append([d_brine.name, '', 1000000., 1000000., mixing_brine_mass])
            if np_round(np_sum(remaining_ion_masses), 3) == 0:
                line = 'DI | ' + format_string + ' g'
            else:
                print('sum of remaining ion masses:', np_sum(remaining_ion_masses))
                line = 'Balance (' + str(np_round(np_sum(remaining_ion_composition), 0)) + ' TDS) | ' + format_string + ' g'
            mixing_text += line.format(total_mass - mixing_brine_mass - added_mass)
            if total_mass - mixing_brine_mass - added_mass != 0.:
                if np_round(np_sum(remaining_ion_masses), 3) == 0:
                    fluids_mixing_list.append(['DI', '', 1000000., 1000000., total_mass - mixing_brine_mass - added_mass])
                else:
                    fluids_mixing_list.append(['Balance (' + str(np_round(np_sum(remaining_ion_composition), 0)) + ' TDS)',
                                               '', 1000000., 1000000., total_mass - mixing_brine_mass - added_mass])
                    self.parent().balance_brine = remaining_ion_composition.tolist()
        formulation_mixing_lines += 2

        # self.final_mixing_label.setText(mixing_text)
        adjusted_mixing_text = ''
        max_len = 0
        for line in mixing_text.split('\n'):
            i = line.index('|')
            if i > max_len:
                max_len = i

        for line in mixing_text.split('\n'):
            i = line.index('|')
            rep_string = ''
            if i < max_len:
                for j in range(max_len - i):
                    rep_string += ' '
            adjusted_mixing_text += line[:i] + rep_string + line[i:]
            adjusted_mixing_text += '\n'

        adjusted_mixing_text = adjusted_mixing_text[:-1]
        self.final_mixing_label.setText(adjusted_mixing_text)
        self.formulation_mixing_lines = formulation_mixing_lines
        self.parent().fluids_mixing_list = fluids_mixing_list
        self.parent().fluids_mixing_salinity = total_fluid_salinity


class Brine(Fluid):

    ion_names = {"Li+": 'lithium', "Na+": 'sodium', "K+": 'potassium', "Mg++": 'magnesium', "Ca++": 'calcium',
                 "Ba++": 'barium', "Sr++": 'strontium', "Fe++": 'ironII', "F-": 'fluoride', "Cl-": 'chloride',
                 "Br-": 'bromide', "HCO3-": 'bicarbonate', "CO3--": 'carbonate', "SO4--": 'sulfate'}

    ion_weights = {"Li+": 6.94, "Na+": 22.99, "K+": 39.1, "Mg++": 24.305, "Ca++": 40.078, "Ba++": 137.327,
                   "Sr++": 87.62, "Fe++": 55.84, "F-": 19.0, "Cl-": 35.453, "Br-": 79.904, "HCO3-": 61.01,
                   "CO3--": 60.008, "SO4--": 96.06}

    additive_salts = {"Na2CO3": [[2, 1], ["Na+", "CO3--"], 105.989], "NaCl": [[1, 1], ["Na+", "Cl-"], 58.443]}
    alkalis = ["Na2CO3"]

    additive_non_salts = {"MEA": 61.08}
    non_salt_alkalis = ["MEA"]

    def __init__(self, name: str, lithium=0., sodium=0., potassium=0., magnesium=0., calcium=0., barium=0.,
                 strontium=0., ironII=0., fluoride=0., chloride=0., bromide=0., bicarbonate=0., carbonate=0.,
                 sulfate=0., sodium_carbonate=0., sodium_chloride=0., mea=0., rv: float=1.):
        super(Brine, self).__init__()
        self.project_manager_list = None
        self.name = name
        self.ref_viscosity = rv  # viscosity at 70 F (21 C)
        self.view_class = BrineTool
        self.composition = {"Li+": lithium, "Na+": sodium, "K+": potassium, "Mg++": magnesium, "Ca++": calcium,
                            "Ba++": barium, "Sr++": strontium, "Fe++": ironII, "F-": fluoride, "Cl-": chloride,
                            "Br-": bromide, "HCO3-": bicarbonate, "CO3--": carbonate, "SO4--": sulfate}
        self.add_salt_composition = {"Na2CO3": sodium_carbonate, "NaCl": sodium_chloride}
        self.add_non_salt_composition = {"MEA": mea}

    @staticmethod
    def get_ordered_keys() -> list:

        return ["Li+", "Na+", "K+", "Mg++", "Ca++", "Ba++", "Sr++", "Fe++",
                "F-", "Cl-", "Br-", "HCO3-", "CO3--", "SO4--"]

    def set_cation(self, cation, value):
        if cation in self.composition.keys():
            self.composition[cation] = value

    def get_tds(self) -> float:

        tds = 0.
        for ppm in self.composition.values():
            tds += ppm
        for percent in self.add_salt_composition.values():
            tds += percent * 10000.

        return tds

    def get_composition_array(self) -> array:

        composition = []
        ordered_keys = Brine.get_ordered_keys()

        for key in ordered_keys:
            composition.append(self.composition[key])

        for key, value in self.add_salt_composition.items():
            add_salt_info = Brine.additive_salts[key]
            for stoic, ion_name in zip(add_salt_info[0], add_salt_info[1]):
                ind = ordered_keys.index(ion_name)
                composition[ind] += value * 10000. * (stoic * Brine.ion_weights[ion_name] / add_salt_info[2])

        return array(composition)

    def get_viscosity(self, temp):
        """From El-Dessouky, Ettouny (2002): Fundamentals of Sea Water Desalination (Appendix A: Themodynamic
        Properties). Variable 'temp' is temperature in degrees Celsius."""

        sal = self.get_tds() / 1000.

        a = 1.474E-3 + 1.5E-5 * temp - 3.927E-8 * power(temp, 2)
        b = 1.0734E-5 - 8.5E-8 * temp + 2.23E-10 * power(temp, 2)
        mur = 1 + a * sal + b * sal ** 2.
        lnmuw = -3.79418 + divide(604.129, 139.18 + temp)

        mu = exp(lnmuw) * mur
        # print('viscosity at ', temp, ' degrees C = ', mu, ' cP')

        return mu

    def get_density(self, temp):
        """From El-Dessouky, Ettouny (2002): Fundamentals of Sea Water Desalination (Appendix A: Themodynamic
        Properties). Variable 'temp' is temperature in degrees Celsius."""

        sal = self.get_tds() / 1000.

        a = (2. * temp - 200.) / 160.
        b = (2. * sal - 150.) / 150.

        f1 = 0.5
        f2 = a
        f3 = 2. * power(a, 2) - 1.
        f4 = 4. * power(a, 2) - 3. * a
        g1 = 0.5
        g2 = b
        g3 = 2. * b ** 2. - 1.

        a1 = 4.032219 * g1 + 0.115313 * g2 + 3.26E-4 * g3
        a2 = -0.108199 * g1 + 1.571E-3 * g2 - 4.23E-4 * g3
        a3 = -0.012247 * g1 + 1.74E-3 * g2 - 9.0E-6 * g3
        a4 = 6.92E-4 * g1 - 8.7E-5 * g2 - 5.3E-5 * g3

        rho = a1 * f1 + a2 * f2 + a3 * f3 + a4 * f4

        return rho

    def get_base_composition(self) -> dict:
        """This method returns a composition for the brine without additive salts."""

        base_composition = self.composition.copy()

        for salt, value in self.add_salt_composition.items():

            salt_ppm = 10000. * value
            ion_data = Brine.additive_salts[salt]
            valences = ion_data[0]
            ions = ion_data[1]
            salt_weight = ion_data[2]

            for ion, valence in zip(ions, valences):

                ion_weight = Brine.ion_weights[ion]
                wt_percent = valence * ion_weight / salt_weight
                ion_ppm = wt_percent * salt_ppm
                base_composition[ion] -= ion_ppm

        return base_composition

    def has_additive_salt(self) -> bool:

        for value in self.add_salt_composition.values():
            if value:
                return True

        return False

    def has_additive_non_salt(self) -> bool:

        for value in self.add_non_salt_composition.values():
            if value:
                return True

        return False

    def n_additive_salts(self) -> int:

        n_add_salt = 0
        for i in range(len(self.add_salt_composition.keys())):
            n_add_salt += 1

        return n_add_salt

    def additive_salts_and_concentrations(self):

        if not self.has_additive_salt():
            return [], []

        salts = []
        concentrations = []

        for salt, c in self.add_salt_composition.items():
            if c > 0.:
                salts.append(salt)
                concentrations.append(c)

        return salts, concentrations

    def additive_non_salts_and_concentrations(self):

        if not self.has_additive_non_salt():
            return [], []

        non_salts = []
        concentrations = []

        for non_salt, c in self.add_non_salt_composition.items():
            if c > 0.:
                non_salts.append(non_salt)
                concentrations.append(c)

        return non_salts, concentrations

    def has_alkali(self) -> bool:

        for alkali in Brine.alkalis:
            if self.add_salt_composition[alkali]:
                return True

        return False

    def n_alkalis(self) -> int:

        n_alkali = 0
        for alkali in Brine.alkalis:
            if self.add_salt_composition[alkali]:
                n_alkali += 1

        return n_alkali

    def first_alkali_and_concentration(self):

        if not self.has_alkali():
            return None, nan

        for alkali in Brine.alkalis:
            if self.add_salt_composition[alkali]:
                return alkali, self.add_salt_composition[alkali]

    def alkalis_and_concentrations(self):

        if not self.has_alkali():
            return [], []

        alkalis = []
        concentrations = []

        for alkali in Brine.alkalis:
            if self.add_salt_composition[alkali]:
                alkalis.append(alkali)
                concentrations.append(self.add_salt_composition[alkali])

        return alkalis, concentrations

    def find_nearest_base_composition_brine(self, brines: list):

        base_composition = self.get_base_composition()
        deviations = []

        for brine in brines:
            total_deviation = 0.
            for ion in base_composition.keys():
                total_deviation += (base_composition[ion] - brine.composition[ion]) ** 2.
            deviations.append(total_deviation)

        min_dev = np_min(deviations)
        i = deviations.index(min_dev)

        return i, brines[i]


class BrineTool(QDialog, utils.SignalView):
    """This class encodes a GUI for entering brine names and compositions into the database brine table,
    and for displaying the names of the brines in the database for selection by the user."""

    def __init__(self, parent, *args):
        """The class constructor creates the DatabaseBrineTool GUI and requires a DatabaseView object
        as input, to be later referenced during information transfer."""
        super(BrineTool, self).__init__(parent=parent)

        # Set simple dialog properties.
        self.setWindowIcon(QIcon('Coreholder Cropped.jpg'))
        self.setWindowTitle('Brine Tool')
        sf = parent.parent().width() / 900.
        self.sf = sf
        self.setFixedSize(int(sf * 470), int(sf * 575))

        self.brine = Brine(name='')  # create a brine object to store composition information
        self.brine_import_data = []
        self.bd_gui = None

        master_layout = QHBoxLayout()   # main layout
        left_widget = QWidget()     # composite widget with name label and edit, composition table, and submit button
        layout = QVBoxLayout()  # layout for left_widget
        left_widget.setLayout(layout)  # set the left_widget layout
        self.setLayout(master_layout)   # set the main layout as the DatabaseBrineTool object's layout

        self.name_label = QLabel(self)      # name label and edit
        self.name_label.setText('Name:')
        self.name_edit = QLineEdit(self)
        layout.addWidget(self.name_label)   # add the name label and edit to the left_widget layout
        layout.addWidget(self.name_edit)

        submit_visible = True
        if args:
            for arg in args:
                if isinstance(arg, Brine):
                    self.brine = arg
                    # print(self.brine)
                if isinstance(arg, str):
                    self.name_edit.setText(arg)
                if isinstance(arg, bool):
                    submit_visible = arg

        table_widget = QWidget(self)
        table_layout = QGridLayout()
        table_widget.setLayout(table_layout)
        self.table = BrineTable(self.brine)   # the brine composition table
        table_layout.addWidget(self.table, 0, 0, 2, 1)    # add the table to the left_widget layout
        self.add_table = BrineTable(self.brine, 'add_salt_composition', 'float')  # the brine add salt composition table
        table_layout.addWidget(self.add_table, 0, 1, 1, 1)  # add the table to the left_widget layout
        self.add_non_table = BrineTable(self.brine, 'add_non_salt_composition', 'float')  # the brine add non-salt composition table
        table_layout.addWidget(self.add_non_table, 1, 1, 1, 1)  # add the table to the left_widget layout
        layout.addWidget(table_widget)

        self.import_button = QPushButton('Import')  # the import button and slot connection
        self.import_button.clicked.connect(self.import_from_db)
        layout.addWidget(self.import_button)  # add the import button to the left_widget layout

        self.submit_button = QPushButton('Submit')  # the submit button and slot connection
        self.submit_button.clicked.connect(self.submit)
        self.submit_button.clicked.connect(self.close)
        layout.addWidget(self.submit_button)    # add the submit button to the left_widget layout
        self.submit_button.setVisible(submit_visible)
        self.name_edit.setReadOnly(not submit_visible)
        self.import_button.setVisible(submit_visible)

        master_layout.addWidget(left_widget)    # add the completed left_widget to the main layout

        # self.brine_list = QListWidget(self)     # make the brine list and connect the slots
        # master_layout.addWidget(self.brine_list)    # add the brine list to the main layout

        self.show()     # show the initialized DatabaseBrineTool

    def import_from_db(self):
        if self.bd_gui is not None:
            return
        from BrineDatabase import BrineDatabaseGUI
        try:
            self.bd_gui = BrineDatabaseGUI(export_app=self)
            self.bd_gui.brine_view.import_button.pack_forget()
            self.bd_gui.brine_view.total_mass_entry.grid_forget()
            self.bd_gui.brine_view.total_mass_label.grid_forget()
            self.bd_gui.brine_view.conc_factor_label.grid_forget()
            self.bd_gui.brine_view.conc_factor_entry.grid_forget()
            self.bd_gui.print_view.grid_forget()
            # bd_width = str(int(self.sf * 680))
            # bd_height = str(int(self.sf * 487))
            # self.bd_gui.geometry(bd_width + 'x' + bd_height)
            self.bd_gui.brine_view.send_button.config(width=20, height=4)
            self.bd_gui.brine_view.button_frame.grid_configure(pady=(15, 16))
            self.bd_gui.brine_view.send_button.config(command=self.export_brine, text='Send to Brine Tool',
                                                      state='normal')
            self.bd_gui.title('Import brine...')
            self.bd_gui.mainloop()
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()
            return

    def export_brine(self):
        brine = self.bd_gui.brine_view.n_brine
        add_salt_composition = brine.add_salt_composition
        brine_from_db = [self.bd_gui.brine_view.abbrev_entry_var.get()]

        for value in brine.composition.values():
            brine_from_db.append(value)

        self.update_brine(brine_from_db, add_salt_composition)

    def update_brine(self, brine_from_db: list, add_salt_composition: dict):
        ion_names_list = ['lithium', 'sodium', 'potassium', 'magnesium', 'calcium', 'barium', 'strontium', 'ironII',
                          'fluoride', 'chloride', 'bromide', 'bicarbonate', 'carbonate', 'sulfate']
        ions = brine_from_db[1:]
        composition = []
        for key in self.brine.composition.keys():
            ind = ion_names_list.index(Brine.ion_names[key])
            composition.append(ions[ind])

        # Fill in the table with table items.
        for i in range(0, self.table.rowCount()):
            self.table.item(i, 0).setText(str(int(composition[i])))

        for i, value in enumerate(add_salt_composition.values()):
            self.add_table.item(i, 0).setText(str(float(value)))

        self.name_edit.setText(brine_from_db[0])

    def submit(self):

        name = self.name_edit.text()

        if not name:
            return

        ref_viscosity = 1.

        composition = []

        r = self.table.rowCount()
        for i in range(r):
            item = self.table.item(i, 0)
            composition.append(float(item.text()))

        r = self.add_table.rowCount()
        for i in range(r):
            item = self.add_table.item(i, 0)
            composition.append(float(item.text()))

        r = self.add_non_table.rowCount()
        for i in range(r):
            item = self.add_non_table.item(i, 0)
            composition.append(float(item.text()))

        self.parent().submit_name(name, *composition, ref_viscosity)
        # self.brine.name = name
        self.brine.ref_viscosity = ref_viscosity

    def clear_ui(self):
        """This function simply clears the DatabaseBrineTool GUI by setting the brine name edit text
        to an empty string and setting all of the compositional values to zeros. The brine list is
        also updated from the database."""

        self.name_edit.setText('')  # set the name edit text to empty string

        for i in range(self.table.rowCount()):
            item = self.table.item(i, 0)    # get item in the table
            if item is not None:
                item.setText('0')   # set the value to 0 if the item exists

    def closeEvent(self, a0: QtGui.QCloseEvent):
        print(self.brine.get_base_composition())
        if self.bd_gui is not None:
            # self.bd_gui.export_app = None
            self.bd_gui.destroy()

        self.close_signal.emit(self)


class BrineTable(QTableWidget):
    """This class augments QTableWidget so that a custom brine table is created, that the user can
    delete brine concentrations using the delete key, and so that input is vetted as positive numeric."""

    def __init__(self, brine: Brine, composition_key: str = 'composition', data_type='int'):
        """The class constructor creates a table of the correct dimensions with correct labels."""
        super(BrineTable, self).__init__()

        self.brine = brine
        self.composition_key = composition_key
        self.data_type = data_type

        # Set the table dimensions.
        composition = brine.__dict__[composition_key].values()
        self.setColumnCount(1)  # only one column
        self.setRowCount(len(composition))

        # Set the column and row labels.
        if composition_key == 'composition':
            self.setHorizontalHeaderLabels(['[ppm]'])
        else:
            self.setHorizontalHeaderLabels(['%'])
        self.setVerticalHeaderLabels(brine.__dict__[composition_key].keys())

        composition = list(composition)
        # Fill in the table with table items.
        for i in range(0, self.rowCount()):
            self.setItem(i, 0, QTableWidgetItem(str(self.convert_to_datatype(composition[i]))))

        self.itemChanged.connect(self.item_changed)

    def convert_to_datatype(self, val):

        if self.data_type == 'float':
            return float(val)
        else:
            return int(val)

    def keyPressEvent(self, e: QKeyEvent):
        """Let the delete key trigger deletion of the ion concentration entry."""
        super(BrineTable, self).keyPressEvent(e)

        if e.key() == Qt.Key_Delete:
            # 'Deleting' is equivalent to zeroing.
            self.currentItem().setText(str(self.convert_to_datatype(0.)))

    def item_changed(self, item: QTableWidgetItem):
        """Make sure new concentrations are positive numeric."""

        try:
            v = float(item.text())
            if v < 0:
                # If the user entered a concentration less than 0, set to 0.
                item.setText(str(self.convert_to_datatype(0.)))
            else:
                self.brine.__dict__[self.composition_key][self.verticalHeaderItem(item.row()).text()] = v
        except ValueError:
            # If the user entered something non-numeric, re-zero.
            item.setText(str(self.convert_to_datatype(0.)))


class Polymer:

    def __init__(self, name: str, pm=None):
        self.name = name
        self.project_manager_list = pm
        self.rheology_solutions = []
        self.view_class = PolymerTool
        self.rheologies_list = utils.SignalListManager()
        self.rheologies_list.add_signal_list()
        self.base_fluids_list = utils.SignalListManager()
        self.base_fluids_list.add_signal_list()

    def get_matching_rheologies(self, brine: Brine) -> list:

        matches = []

        for pr in self.rheologies_list.signal_lists[0].objects:
            if isinstance(pr.base_fluid.specific_fluid.brine, Brine):
                if pr.base_fluid.specific_fluid.brine == brine:
                    matches.append(pr)

        return matches

    def get_viscosity_vs_concentration_data_interp(self, brine: Brine, temp: float, shear: float):

        matches = self.get_matching_rheologies(brine)
        if not matches:
            return [], []

        for match in matches:
            cs, vs = match.get_viscosity_vs_concentration_data_interp(temp, shear)
            if not any(isnan(vs)):
                return cs, vs

        return [], []

    def estimate_viscosity(self, brine: Brine, temp: float, shear: float, conc: float, show: bool=False):

        matches = self.get_matching_rheologies(brine)

        if not matches:
            return nan

        rheo_to_use = matches[-1]
        return rheo_to_use.estimate_viscosity(temp, shear, conc, show)

    def estimate_concentration(self, viscosity: float, brine: Brine, temp: float, shear: float):

        matches = []

        for r in self.rheology_solutions:
            assert isinstance(r, InjectionFluid)
            sf = r.specific_fluid
            assert isinstance(sf, BrineInjectionFluid)
            if sf.brine == brine:
                matches.append(r)

        if not matches:
            return nan

        concentrations = []
        viscosities = []

        for match in matches:
            ps = match.specific_fluid.polymer_solution
            concentrations.append(ps.concentration)
            # print(ps, len(ps.stored_averages))
            viscosities.append(ps.get_viscosity(temp, shear))

        if len(concentrations) == 1:
            return concentrations[0]
        elif len(concentrations) == 2:
            m = (concentrations[1] - concentrations[0]) / (viscosities[1] - viscosities[0])
            return viscosities[0] + m * (viscosity - viscosities[0])
        elif len(concentrations) == 3:
            x1 = concentrations[0]
            x2 = concentrations[1]
            x3 = concentrations[2]
            y1 = viscosities[0]
            y2 = viscosities[1]
            y3 = viscosities[2]

            dx21 = x2 - x1
            dx31 = x3 - x1
            dx221 = x2 ** 2. - x1 ** 2.
            dx231 = x3 ** 2. - x1 ** 2.
            dy21 = y2 - y1
            dy31 = y3 - y1

            a = (dy31 * dx21 - dy21 * dx31) / (dx21 * dx231 - dx221 * dx31)
            b = (dy21 - a * dx221) / dx21
            c = y1 - a * x1 ** 2. - b * x1

            return (-b + sqrt(b ** 2. - 4 * a * c)) / (2. * a)

        # print(concentrations, viscosities)
        return 1.


class PolymerTool(QDialog, utils.SignalView):

    def __init__(self, poly: Polymer, parent, *args):
        # print('Polymer Tool:', poly, parent, args)
        super(PolymerTool, self).__init__(parent=parent)
        sf = parent.width() / 900.
        self.sf = sf
        self.setFixedSize(int(sf * 600), int(sf * 300))
        self.polymer = poly
        self.setWindowTitle(poly.name + ' Viewer')
        self.polymer_rheology_view = None

        lyt = QGridLayout()
        self.setLayout(lyt)

        list_name = ['Base Fluids']
        cls = [[BrineInjectionFluidRheologyView], [InjectionFluid]]
        args = [None]
        self.slm_bf = self.polymer.base_fluids_list
        self.list_widget_bf = utils.SignalListManagerWidget(self, self.polymer.base_fluids_list, list_name, cls, *args)
        self.slm_bf.add_signal_list()
        for bf in poly.base_fluids_list.signal_lists[0].objects:
            bf.project_manager_list = self.list_widget_bf

        list_name = ['Rheologies']
        cls = [[polymer_rheology_view_wrapper], [PolymerRheology]]
        args = [['Name'], [poly]]
        self.slm = self.polymer.rheologies_list
        self.list_widget = utils.SignalListManagerWidget(self, self.polymer.rheologies_list, list_name, cls, *args)
        self.slm.add_signal_list()
        for rheo in poly.rheologies_list.signal_lists[0].objects:
            rheo.project_manager_list = self.list_widget

        pmv = self.parent()

        self.list_widget_bf.pm_list.signal_lists[0].set_ref_lists([pmv.brine_list.signal_lists[0]])
        self.list_widget.pm_list.signal_lists[0].set_ref_lists([self.list_widget_bf.pm_list.signal_lists[0]])
        # self.slm.signal_lists[0].set_ref_lists([self.slm_bf.signal_lists[0]])

        lyt.addWidget(self.list_widget_bf, 0, 0, 1, 1)
        lyt.addWidget(self.list_widget, 0, 1, 1, 1)

        self.show()

    def closeEvent(self, a0: QtGui.QCloseEvent):
        if self.polymer_rheology_view is not None:
            self.polymer_rheology_view.close()
        self.close_signal.emit(self)


class PolymerRheology:

    def __init__(self, name: str, ref_objects: list, *args):
        self.name = name
        self.ref_objects = ref_objects
        self.polymer_solutions_list = []

        everything = [*ref_objects, *args]
        self.poly = None
        for thing in everything:
            if isinstance(thing, Polymer):
                self.poly = thing
        if ref_objects:
            self.base_fluid = ref_objects[0]
        else:
            self.base_fluid = None

        self.view_class = PolymerRheologyView

    def add_rheology(self, parent):
        try:
            ps = PolymerSolution(base_fluid=self.base_fluid.specific_fluid.brine, primary_additive=self.poly)
            tool = PolymerRheologyImportTool(ps, parent, self)
            return tool
        except Exception as e:
            print(e)

    def get_viscosity_vs_concentration_data_interp(self, temp: float, shear: float):

        cs = []
        vs = []

        for ps in self.polymer_solutions_list:
            cs.append(ps.concentration)
            vs.append(ps.get_viscosity_interp(temp, shear))

        return cs, vs

    def estimate_viscosity(self, temp: float, shear: float, conc: float, show: bool=False):

        cs = []
        ts = []
        ss = []
        ctss = []
        vs = []

        for ps in self.polymer_solutions_list:
            for rt, sm, sa, sd in zip(ps.stored_ref_temps, ps.stored_models, ps.stored_averages, ps.stored_shear_data):
                for s in sd:
                    cs.append(ps.concentration)
                    ts.append(rt)
                    ss.append(s)
                    ctss.append([ps.concentration, rt, s])
                    ps.model = deepcopy(sm)
                    ps.ref_temp = rt
                    if not isnan(sm['n']):
                        vs.append(ps.rheology_model(s, sm['eta_o'], sm['eta_inf'], sm['n'], sm['l']))
                    else:
                        vs.append(sa)

        if not vs:
            return nan

        if len(cs) == 1 and cs[0] != conc:
            return nan

        cs = array(cs)
        ss = array(ss)
        vs = array(vs)
        ts = array(ts)
        ctss = array(ctss)

        try:
            params = [-0.01274, 0.807189, 20.87, 0.049751, 0.025543, 0.020635]
            popt, _ = curve_fit(self.jouenne_model_bundle, ctss, vs, p0=params,
                                bounds=([-0.1, 0., 0.1, 0., 0., 0.], [0., 2., 100., 1., 1., 1.]))
            v_est = self.jouenne_model(conc, temp, shear, *popt)

            if show:

                fig = plt.figure()
                ax = fig.add_subplot(projection='3d')
                ax.scatter(cs, log10(ss), log10(vs), c=ts, cmap='magma', vmin=20., vmax=100.)
                cmap = cm.get_cmap('magma')
                rgba = cmap((temp - 20.) / 80.)
                rgba = (*rgba[:3], 0.3)

                c_lin = linspace(floor(np_min(cs)), ceil(np_max(cs)), 40)
                s_lin = linspace(-1, 3, 40)
                c_mesh, s_mesh = meshgrid(c_lin, s_lin)
                t_mesh = temp * ones(shape=c_mesh.shape)
                v_mesh = self.jouenne_model(c_mesh, t_mesh, power(10., s_mesh), *popt)
                # v_mesh_est = v_est * ones(shape=c_mesh.shape)

                ax.plot_wireframe(c_mesh, s_mesh, log10(v_mesh), color=rgba)
                # ax.plot_wireframe(c_mesh, s_mesh, log10(v_mesh_est), color=(0., 0., 0., 0.3))
                ax.plot(array([conc, conc]), log10(array([shear, shear])), array([0., 3.]), color='red')
                ax.plot(array([conc, conc]), log10(array([0.1, 1000.])), log10(array([v_est, v_est])),
                        color='red')
                ax.plot(array([floor(np_min(cs)), ceil(np_max(cs))]), log10(array([shear, shear])),
                        log10(array([v_est, v_est])), color='red')

                plt.xlabel('[ppm]')
                plt.ylabel(r'$log_{10} [1/s]$')
                title_str = 'V_est = %.1f cP' % v_est
                plt.title(title_str)

                plt.show()

        except Exception as e:
            v_est = nan
            print(e)

        return v_est

    def jouenne_model(self, cs, ts, ss, *params):

        n_m = params[0]
        n_b = params[1]
        iv = params[2]
        iv_b1 = params[3]
        iv_b2 = params[4]
        l_d = params[5]

        bv = self.base_fluid.specific_fluid.brine.get_viscosity(ts)
        bd = self.base_fluid.specific_fluid.brine.get_density(ts)

        c = 1.e-4 * divide(cs, bd) # [g/dL]
        eta0 = multiply(bv, 1. + c * iv + iv_b1 * power(c * iv, 2) + iv_b2 * power(c * iv, 3))
        l = l_d * divide(np_min(ts) + 273.15, ts + 273.15) * (1. + 0.04 * power(c * iv, 2.4))
        n = n_b + n_m * c * iv

        vs = bv + divide(eta0 - bv, power(1. + multiply(l, ss), 0.5 - 0.5 * n))

        return vs

    def jouenne_model_bundle(self, ctss, *params):

        cs = ctss[:, 0]
        ts = ctss[:, 1]
        ss = ctss[:, 2]

        return self.jouenne_model(cs, ts, ss, *params)


def polymer_rheology_view_wrapper(parent, *args):
    PolymerRheologyView(None, parent, *args)


class PolymerRheologyEdit(QLineEdit):

    def __init__(self, parent, min_val: float=0.):
        super(PolymerRheologyEdit, self).__init__(parent=parent)
        self.min_val = min_val
        self.current_val = nan

        self.editingFinished.connect(self.check_val)
        self.editingFinished.connect(parent.calculate_viscosity)

    def check_val(self):

        txt = self.text()
        try:
            val = float(txt)
            if val < self.min_val:
                raise ValueError
            self.current_val = val
        except ValueError:
            if not isnan(self.current_val):
                self.setText(str(self.current_val))
            else:
                self.setText('')


class PolymerRheologyView(QDialog):

    def __init__(self, p_rheo, parent, *args):
        super(PolymerRheologyView, self).__init__(parent=parent)
        parent.polymer_rheology_view = self
        self.import_tool = None
        # print('PolymerRheologyView args:', args)

        if isinstance(parent.parent(), PolymerTool):
            sf = parent.parent().width() / 600.
            self.pmv = parent.parent().parent()
            # print('parent.parent() is PolymerTool', self.pmv)
            sl = parent.parent().slm_bf.signal_lists[0]
            cr = sl.current_row
            if cr == -1:
                return
            ref_objs = [sl.objects[cr]]
        else:
            ref_objs = []
            sf = parent.parent().width() / 900.
            self.pmv = parent.parent()
            # print('parent.parent() is not PolymerTool', self.pmv)

        self.setFixedSize(int(sf * 200), int(sf * 100))
        lyt = QGridLayout()
        self.setLayout(lyt)
        self.name_edit = QLineEdit(parent=self)

        if p_rheo is None:
            self.polymer_rheology = PolymerRheology(args[0][0], ref_objs, args[1][0])
            self.submit_button = QPushButton(parent=self, text='Submit')
            self.submit_button.clicked.connect(self.submit_button_clicked)
            lyt.addWidget(self.name_edit, 0, 0, 1, 1)
            lyt.addWidget(self.submit_button, 1, 0, 1, 1)
            self.t_edit = None
            self.c_edit = None
            self.s_edit = None
        else:
            self.polymer_rheology = p_rheo
            self.name_edit.setText(p_rheo.name)
            self.import_button = QPushButton(parent=self, text='Import Rheology')
            self.import_button.clicked.connect(self.import_button_clicked)
            t_label = QLabel(parent=self, text='T [C]: ')
            self.t_edit = PolymerRheologyEdit(parent=self)
            c_label = QLabel(parent=self, text='c [ppm]: ')
            self.c_edit = PolymerRheologyEdit(parent=self)
            s_label = QLabel(parent=self, text=chr(947) + u' [1/s]: ')
            self.s_edit = PolymerRheologyEdit(parent=self, min_val=0.1)
            self.viscosity_label = QLabel(parent=self)
            self.view_button = QPushButton(parent=self, text='View')
            self.view_button.clicked.connect(self.view_button_clicked)
            if len(self.polymer_rheology.polymer_solutions_list) < 2:
                self.view_button.setEnabled(False)
            self.setFixedSize(int(sf * 200), int(sf * 175))
            lyt.addWidget(self.name_edit, 0, 0, 1, 4)
            lyt.addWidget(self.import_button, 1, 0, 1, 4)
            lyt.addWidget(t_label, 2, 0, 1, 1)
            lyt.addWidget(self.t_edit, 2, 1, 1, 1)
            lyt.addWidget(c_label, 3, 0, 1, 1)
            lyt.addWidget(self.c_edit, 3, 1, 1, 1)
            lyt.addWidget(s_label, 4, 0, 1, 1)
            lyt.addWidget(self.s_edit, 4, 1, 1, 1)
            lyt.addWidget(self.viscosity_label, 2, 2, 1, 2)
            lyt.addWidget(self.view_button, 3, 2, 2, 2)

            self.refresh_viscosity_estimation()

        if ref_objs:
            self.polymer_rheology.base_fluid = ref_objs[0]

        self.show()

    def refresh_viscosity_estimation(self):
        if len(self.polymer_rheology.polymer_solutions_list) < 2:
            self.enable_estimation(False)
        else:
            self.enable_estimation(True)
        self.viscosity_label.setText('')

    def enable_estimation(self, enable: bool):

        self.t_edit.setEnabled(enable)
        self.c_edit.setEnabled(enable)
        self.s_edit.setEnabled(enable)

    def calculate_viscosity(self, plot_me: bool=False):

        pr = self.polymer_rheology
        if len(pr.polymer_solutions_list) < 2:
            self.viscosity_label.setText('')
            return

        temp = self.t_edit.text()
        c = self.c_edit.text()
        s = self.s_edit.text()

        if not temp or not c or not s:
            return
        else:
            temp = float(temp)
            c = float(c)
            s = float(s)

        try:
            v_est = pr.poly.estimate_viscosity(pr.base_fluid.specific_fluid.brine, temp, s, c, plot_me)
            v_text = '%.1f cP' % v_est
            self.viscosity_label.setText(v_text)
            cs, vs = pr.get_viscosity_vs_concentration_data_interp(temp, s)
            print('concentration data:', cs, vs)
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()

    def submit_button_clicked(self):

        name = self.name_edit.text()
        if not name:
            return

        if self.polymer_rheology.base_fluid is None:
            return

        bf = self.polymer_rheology.base_fluid
        args = [self.polymer_rheology.poly, bf]

        # print(name, bf, *args)

        # print(self.parent().parent())
        self.parent().submit_name_ref_objects(name, [bf], *args)
        self.close()

    def import_button_clicked(self):

        if self.import_tool is None:
            self.import_tool = self.polymer_rheology.add_rheology(self)

    def view_button_clicked(self):

        self.calculate_viscosity(True)

    def closeEvent(self, a0: QtGui.QCloseEvent):

        if self.import_tool is not None:
            self.import_tool.close()

        if self.t_edit is not None:
            self.t_edit.blockSignals(True)
            self.c_edit.blockSignals(True)
            self.s_edit.blockSignals(True)

        try:
            self.parent().polymer_rheology_view = None
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()


class BrineInjectionFluidRheologyView(SpecificFluidView):

    def __init__(self, parent, brine_if: InjectionFluid=None):
        super(BrineInjectionFluidRheologyView, self).__init__(parent=parent, injection_fluid=brine_if)

        # this is probably built on an incorrect premise...the "else" likely executes all the time, and works in the
        # subsequent code
        # if isinstance(self.parent(), PolymerTool):
        #     self.pml = self.parent().parent().parent()
        # else:
        #     self.pml = self.parent().parent().parent().parent()

        if isinstance(self.parent().parent(), PolymerTool):
            self.pml = self.parent().parent().parent()
        else:
            self.pml = self.parent().parent().parent().parent()

        sf = self.pml.width() / 900.
        self.sf = sf
        # self.setFixedWidth(int(sf * 300))
        self.setFixedSize(int(sf * 300), int(sf * 100))
        self.add_button.setFixedSize(int(sf * 25), int(sf * 25))
        self.remove_button.setFixedSize(int(sf * 25), int(sf * 25))
        self.submit_button.setFixedHeight(int(sf * 25))
        self.name_text.setText('Base Fluid Name:')
        source_names = self.pml.brine_list.signal_lists[0].item_names
        self.source_fluid_combo_box.addItems(source_names)

        mode = 'build'
        if brine_if is not None:
            mode = 'view'

        self.setWindowTitle('Base Fluid for Polymer Rheology')

        if mode == 'view':

            brine_objects = self.pml.brine_list.signal_lists[0].objects
            for i, brine_object in enumerate(brine_objects):
                if brine_object == brine_if.ref_objects[0]:
                    self.source_fluid_combo_box.setCurrentIndex(i)
            if len(brine_if.ref_objects) > 1:
                additives = brine_if.ref_objects[1:]
                for i, additive in enumerate(additives):

                    self.add_additive(i)
                    w = self.combo_widgets[-1]
                    w.blockSignals(True)

                    index = 0
                    if isinstance(additive, CoSolvent):
                        index = 0
                    elif isinstance(additive, Surfactant):
                        index = 1
                    elif isinstance(additive, Scavenger):
                        index = 2
                    elif isinstance(additive, Other):
                        index = 3

                    w.combo_box.setCurrentIndex(index)
                    additive_stocks = w.source_lists[0].signal_lists[index + 1].objects
                    for j, additive_stock in enumerate(additive_stocks):
                        if additive_stock == additive:
                            w.combo_box_2.setCurrentIndex(j+1)
                            w.edit.setText(str(int(brine_if.concentrations[i])))

                    w.blockSignals(False)

            # if isinstance(brine_if.specific_fluid, BrineInjectionFluid):
            self.fluid_name_edit.setText(brine_if.name)
            #     if brine_if.specific_fluid.polymer_solution is not None:
            #         if isinstance(self.parent().parent(), PolymerTool):
            #             PolymerRheologyImportTool(brine_if.specific_fluid.polymer_solution, self.parent())
            #         else:
            #             PolymerSolutionView(brine_if.specific_fluid.polymer_solution, self.parent())

        self.show()

    def add_additive(self, conc_index: int=-1):
        new_widget = AdditiveRheologyWidget([self.pml.chemical_list], self, conc_index)
        self.layout().addWidget(new_widget)
        if self.injection_fluid is not None:
            new_widget.edit.editingFinished.connect(new_widget.update_injection_fluid_concentration)
            new_widget.combo_box.setEnabled(False)
            new_widget.combo_box_2.setEnabled(False)
        self.combo_widgets.append(new_widget)
        self.setFixedHeight(self.height() + int(self.sf * 51))

    def submit(self):

        name = self.fluid_name_edit.text()
        if not name:
            return

        k = self.source_fluid_combo_box.currentIndex()
        brine_obj = self.pml.brine_list.signal_lists[0].objects[k]

        concentrations = []
        ref_objs = [brine_obj]
        for w in self.combo_widgets:
            i = w.combo_box.currentIndex()
            j = w.combo_box_2.currentIndex()
            if j > 0:
                ref_obj = w.source_lists[0].signal_lists[i+1].objects[j-1]
                ref_objs.append(ref_obj)
            concentrations.append(float(w.edit.text()))

        args = [[concentrations, 'brine-based']]

        try:
            success = self.parent().submit_name_ref_objects(name, ref_objs, *args)
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()
            try:
                self.parent().submit_name_ref_objects(name, ref_objs, *args)
            except Exception as ex:
                QMessageBox(parent=self, text=str(ex)).exec_()
            return

        if success:
            bf = self.parent().polymer_rheology.base_fluid
            try:
                BrineInjectionFluidRheologyView(self.parent(), bf)
                self.close()
            except Exception as e:
                QMessageBox(parent=self, text=str(e)).exec_()


class AdditiveRheologyWidget(QWidget):

    def __init__(self, source_lists: list, parent=None, conc_index=-1):
        super(AdditiveRheologyWidget, self).__init__(parent=parent)

        self.source_lists = source_lists
        self.conc_index = conc_index

        h_layout = QHBoxLayout()
        h_layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(h_layout)
        self.combo_box = QComboBox(parent=self)
        self.combo_box.addItems(['Co-solvent', 'Surfactant', 'Scavenger', 'Other'])
        self.combo_box.currentIndexChanged.connect(self.set_combobox_names)
        self.combo_box_2 = QComboBox(parent=self)

        self.edit = QLineEdit(parent=self)
        self.edit.setFixedWidth(int(parent.sf * 60))
        text = QLabel(parent=self, text='[ppm]')
        text.setFixedWidth(int(parent.sf * 40))
        h_layout.addWidget(self.combo_box)
        h_layout.addWidget(self.combo_box_2)
        h_layout.addWidget(self.edit)
        h_layout.addWidget(text)

        self.set_combobox_names()

    def set_combobox_names(self):
        self.combo_box_2.clear()
        cb_list = list(' ')
        i = self.combo_box.currentIndex()
        names = self.source_lists[0].signal_lists[i + 1].item_names
        for i, name in enumerate(names):
            cb_list.append(name)
        self.combo_box_2.addItems(cb_list)

    def update_injection_fluid_concentration(self):

        if self.parent().injection_fluid is None:
            return

        if self.conc_index == -1:
            return

        # print(self.conc_index, float(self.edit.text()))
        # print('concentrations before: ', self.parent().injection_fluid.concentrations)
        self.parent().injection_fluid.concentrations[self.conc_index] = float(self.edit.text())
        # if self.combo_box.currentIndex() == 0 and self.is_brine:
        #     self.parent().injection_fluid.specific_fluid.polymer_solution.concentration = float(self.edit.text())
        # print('concentrations after: ', self.parent().injection_fluid.concentrations)


class CoSolvent:

    def __init__(self, name: str, pm=None):
        self.name = name
        self.project_manager_list = pm


class Surfactant:

    def __init__(self, name: str, pm=None):
        self.name = name
        self.project_manager_list = pm


class Formulation:

    def __init__(self, name: str, ref_objects: list, *args):
        print('Formulation __init__', name, ref_objects, args)
        # super(Formulation, self).__init__(name=name, ref_objects=ref_objects, *args)

        self.brine_if = InjectionFluid(name, ref_objects, *args)
        print('Formulation brine_if:', self.brine_if)

        refs_valid = True
        if not len(ref_objects) > 1:
            refs_valid = False

        has_surfactant = False
        for ref in ref_objects[1:]:
            if isinstance(ref, Surfactant):
                has_surfactant = True
                break

        if not has_surfactant:
            refs_valid = False

        if not refs_valid:
            raise ValueError('Formulation ref_objects must contain a brine and at least one surfactant.')

        self.name = name
        self.project_manager_list = None
        self.ref_objects = ref_objects
        self.view_class = FormulationView
        self.sc_list = ref_objects[1:]
        self.concentrations = args[0][0]

        print(self.name, self.brine_if, self.sc_list, self.concentrations)

    def has_cosolvent(self) -> bool:

        for item in self.sc_list:
            if isinstance(item, CoSolvent):
                return True

        return False


def formulation_view(formulation: Formulation, parent):
    BrineInjectionFluidFormulationView(parent=parent, brine_if=formulation.brine_if)


class FormulationView(QObject, utils.SignalView):

    def __init__(self, formulation: Formulation, parent):
        super(FormulationView, self).__init__()
        real_view = BrineInjectionFluidFormulationView(parent=parent, brine_if=formulation.brine_if)
        real_view.close_signal.connect(self.relay_close_signal)

    def relay_close_signal(self, _):
        self.close_signal.emit(self)


class BrineInjectionFluidFormulationViewWrapper:

    def __init__(self, parent, *args):
        BrineInjectionFluidFormulationView(parent)


class BrineInjectionFluidFormulationView(SpecificFluidView):

    def __init__(self, parent, brine_if: InjectionFluid=None):
        super(BrineInjectionFluidFormulationView, self).__init__(parent=parent, injection_fluid=brine_if)

        # this is probably built on an incorrect premise...the "else" likely executes all the time, and works in the
        # subsequent code
        # if isinstance(self.parent(), PolymerTool):
        #     self.pml = self.parent().parent().parent()
        # else:
        #     self.pml = self.parent().parent().parent().parent()

        if isinstance(parent, utils.SignalListManagerWidget):
            self.pml = parent.parent()
        else:
            print(parent)
            self.pml = parent

        sf = self.pml.width() / 900.
        self.sf = sf
        # self.setFixedWidth(int(sf * 300))
        self.setFixedSize(int(sf * 300), int(sf * 100))
        self.add_button.setFixedSize(int(sf * 25), int(sf * 25))
        self.remove_button.setFixedSize(int(sf * 25), int(sf * 25))
        self.submit_button.setFixedHeight(int(sf * 25))
        self.name_text.setText('Formulation Name:')
        source_names = self.pml.brine_list.signal_lists[0].item_names
        self.source_fluid_combo_box.addItems(source_names)

        mode = 'build'
        if brine_if is not None:
            mode = 'view'

        self.setWindowTitle('Formulation Composition')

        if mode == 'view':

            brine_objects = self.pml.brine_list.signal_lists[0].objects
            for i, brine_object in enumerate(brine_objects):
                if brine_object == brine_if.ref_objects[0]:
                    self.source_fluid_combo_box.setCurrentIndex(i)
            if len(brine_if.ref_objects) > 1:
                additives = brine_if.ref_objects[1:]
                for i, additive in enumerate(additives):

                    self.add_additive(i)
                    w = self.combo_widgets[-1]
                    w.blockSignals(True)

                    index = 0
                    if isinstance(additive, CoSolvent):
                        index = 0
                    elif isinstance(additive, Surfactant):
                        index = 1

                    w.combo_box.setCurrentIndex(index)
                    additive_stocks = w.source_lists[0].signal_lists[index + 1].objects
                    for j, additive_stock in enumerate(additive_stocks):
                        if additive_stock == additive:
                            w.combo_box_2.setCurrentIndex(j+1)
                            w.edit.setText(str(int(brine_if.concentrations[i])))

                    w.blockSignals(False)

            # if isinstance(brine_if.specific_fluid, BrineInjectionFluid):
            self.fluid_name_edit.setText(brine_if.name)
            #     if brine_if.specific_fluid.polymer_solution is not None:
            #         if isinstance(self.parent().parent(), PolymerTool):
            #             PolymerRheologyImportTool(brine_if.specific_fluid.polymer_solution, self.parent())
            #         else:
            #             PolymerSolutionView(brine_if.specific_fluid.polymer_solution, self.parent())

        self.show()

    def add_additive(self, conc_index: int=-1):
        new_widget = AdditiveFormulationWidget([self.pml.chemical_list], self, conc_index)
        self.layout().addWidget(new_widget)
        if self.injection_fluid is not None:
            new_widget.edit.editingFinished.connect(new_widget.update_injection_fluid_concentration)
            new_widget.combo_box.setEnabled(False)
            new_widget.combo_box_2.setEnabled(False)
        self.combo_widgets.append(new_widget)
        self.setFixedHeight(self.height() + int(self.sf * 51))

    def submit(self):

        name = self.fluid_name_edit.text()
        if not name:
            return

        k = self.source_fluid_combo_box.currentIndex()
        brine_obj = self.pml.brine_list.signal_lists[0].objects[k]

        concentrations = []
        ref_objs = [brine_obj]
        for w in self.combo_widgets:
            i = w.combo_box.currentIndex()
            j = w.combo_box_2.currentIndex()
            if j > 0:
                ref_obj = w.source_lists[0].signal_lists[i+1].objects[j-1]
                ref_objs.append(ref_obj)
            concentrations.append(float(w.edit.text()))

        args = [[concentrations, 'brine-based']]

        try:
            self.parent().submit_name_ref_objects(name, ref_objs, *args)
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()
            try:
                self.parent().submit_name_ref_objects(name, ref_objs, *args)
            except Exception as ex:
                QMessageBox(parent=self, text=str(ex)).exec_()
            return

        # if success:
        #     bf = self.parent().polymer_rheology.base_fluid
        #     try:
        #         BrineInjectionFluidFormulationView(self.parent(), bf)
        #         self.close()
        #     except Exception as e:
        #         QMessageBox(parent=self, text=str(e)).exec_()


class AdditiveFormulationWidget(QWidget):

    def __init__(self, source_lists: list, parent=None, conc_index=-1):
        super(AdditiveFormulationWidget, self).__init__(parent=parent)

        self.source_lists = source_lists
        self.conc_index = conc_index

        h_layout = QHBoxLayout()
        h_layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(h_layout)
        self.combo_box = QComboBox(parent=self)
        self.combo_box.addItems(['Co-solvent', 'Surfactant'])
        self.combo_box.currentIndexChanged.connect(self.set_combobox_names)
        self.combo_box_2 = QComboBox(parent=self)

        self.edit = QLineEdit(parent=self)
        self.edit.setFixedWidth(int(parent.sf * 60))
        text = QLabel(parent=self, text='[ppm]')
        text.setFixedWidth(int(parent.sf * 40))
        h_layout.addWidget(self.combo_box)
        h_layout.addWidget(self.combo_box_2)
        h_layout.addWidget(self.edit)
        h_layout.addWidget(text)

        self.set_combobox_names()

    def set_combobox_names(self):
        self.combo_box_2.clear()
        cb_list = list(' ')
        i = self.combo_box.currentIndex()
        names = self.source_lists[0].signal_lists[i + 1].item_names
        for i, name in enumerate(names):
            cb_list.append(name)
        self.combo_box_2.addItems(cb_list)

    def update_injection_fluid_concentration(self):

        if self.parent().injection_fluid is None:
            return

        if self.conc_index == -1:
            return

        # print(self.conc_index, float(self.edit.text()))
        # print('concentrations before: ', self.parent().injection_fluid.concentrations)
        self.parent().injection_fluid.concentrations[self.conc_index] = float(self.edit.text())
        # if self.combo_box.currentIndex() == 0 and self.is_brine:
        #     self.parent().injection_fluid.specific_fluid.polymer_solution.concentration = float(self.edit.text())
        # print('concentrations after: ', self.parent().injection_fluid.concentrations)


class Scavenger:

    def __init__(self, name: str, pm=None):
        self.name = name
        self.project_manager_list = pm


class Other:

    def __init__(self, name: str, pm=None):
        self.name = name
        self.project_manager_list = pm


class PolymerSolution(Fluid):

    def __init__(self, base_fluid, primary_additive, concentration=0):
        super(PolymerSolution, self).__init__()

        self.base_fluid = base_fluid
        self.primary_additive = primary_additive

        self.averaging_on = False
        self.average_viscosity = nan
        self.stored_averages = []

        self.concentration = concentration
        self.model = {"eta_o": nan, "eta_inf": nan, "n": nan, "l": nan}
        self.stored_models = []
        self.ref_temp = nan
        self.stored_ref_temps = []
        self.shear_data = nan
        self.stored_shear_data = []
        self.viscosity_data = nan
        self.stored_viscosity_data = []
        self.test_name = ''
        self.stored_test_names = []
        self.excluded = []
        self.stored_excluded = []

    def clear_stored(self):
        self.stored_averages = []
        self.stored_models = []
        self.stored_ref_temps = []
        self.stored_shear_data = []
        self.stored_viscosity_data = []
        self.stored_test_names = []
        self.stored_excluded = []

    def initialize_stored(self, i: int, cl: bool=True):

        if cl:
            self.clear_stored()

        for j in range(i):
            self.stored_averages.append(nan)
            self.stored_models.append({'eta_o': nan, 'eta_inf': nan, 'n': nan, 'l': nan})
            self.stored_ref_temps.append(nan)
            self.stored_shear_data.append([])
            self.stored_viscosity_data.append([])
            self.stored_test_names.append([])
            self.stored_excluded.append([])

    def get_viscosity_interp(self, temp: float, shear: float) -> float:

        temp = np_round(temp, 0)
        ref_temps = np_round(array(self.stored_ref_temps), 0)
        ref_temps = list(ref_temps)

        if temp not in ref_temps:
            return nan

        index = ref_temps.index(temp)
        sd = array(self.stored_shear_data[index])
        order = sd.argsort()
        log_sd = log10(sd)
        vd = array(self.stored_viscosity_data[index])
        log_vd = log10(vd)

        return power(10., interp(log10(shear), log_sd[order], log_vd[order]))

    def get_viscosity(self, temp, shear=0.):
        # this function returns the shear viscosity of the polymer solution at the current
        # temperature and shear rate
        # print(self.base_fluid)

        if isinstance(self.base_fluid, Brine):
            num_no_nan = where(~isnan(self.stored_averages))[0].size
            print(self.stored_averages, self.stored_ref_temps, self.stored_models)
            if num_no_nan == 1:
                print('one viscosity measurement for polymer.')
                brine_viscosity = self.base_fluid.get_viscosity(temp)
                brine_viscosity_ref = self.base_fluid.get_viscosity(self.stored_ref_temps[0])
                self.model = deepcopy(self.stored_models[0])
                eta_p = self.get_reference_viscosity(shear)
                # print('Polymer relative viscosity: ', eta_p, ' cP')

                return (brine_viscosity / brine_viscosity_ref) * eta_p  # viscosity [cP]
            else:
                print('more than one viscosity measurement for polymer.')
                # brine_viscosity = self.base_fluid.get_viscosity(temp)
                # brine_viscosity_ref = self.base_fluid.get_viscosity((9. / 5.) * self.ref_temp + 32.)
                # eta_p = self.get_reference_viscosity(shear)

                if self.averaging_on:
                    plt.plot(self.stored_ref_temps, self.stored_averages)
                    plt.show()
                    temps = 9 * array(self.stored_ref_temps) / 5 + 32
                    viscosities = array(self.stored_averages)
                    fit_params, _ = curve_fit(PolymerSolution.temp_dependence_polymer, temps, viscosities,
                                              bounds=(0, 10000))
                    temp_o = fit_params[0]
                    eta = PolymerSolution.temp_dependence_polymer(temp, temp_o)
                    # print(round(eta, 2), ' [cP]')
                    return eta
                else:
                    viscosities = []
                    for ref_temp, model in zip(self.stored_ref_temps, self.stored_models):
                        # print(ref_temp, model)
                        visc = PolymerSolution.rheology_model(shear,
                                                              model['eta_o'], model['eta_inf'], model['n'], model['l'])
                        viscosities.append(visc)

                    temps = array(self.stored_ref_temps)
                    viscosities = array(viscosities)
                    # plt.plot(temps, viscosities)
                    # plt.show()
                    fit_params, _ = curve_fit(PolymerSolution.temp_dependence_polymer, temps, viscosities,
                                              bounds=(0, 10000))
                    temp_o = fit_params[0]
                    eta = PolymerSolution.temp_dependence_polymer(temp, temp_o)
                    # print(round(eta, 2), ' [cP]')
                    return eta

        elif isinstance(self.base_fluid, OilSample):
            # print('calculating viscosity for OilSample')
            # print('ref_temps: ', self.stored_ref_temps)
            # print('averages: ', self.stored_averages)
            temps = array(self.stored_ref_temps)
            num_no_nan = where(~isnan(temps))[0].size
            viscosities = array(self.stored_averages)

            if num_no_nan > 1:
                fit_params = linregress(temps, log10(log10(viscosities)))
                # print(fit_params)
                m = fit_params[0]
                b = fit_params[1]

            else:
                m = 0.
                b = log10(log10(viscosities[0]))

            eta = PolymerSolution.temp_dependence_oil(temp, m, b)

            # print(round(eta, 2), ' [cP]')
            return eta

            # return (brine_viscosity / brine_viscosity_ref) * eta_p  # viscosity [cP]

    def get_reference_viscosity(self, shear=0.):
        # this function returns the viscosity at the reference temperature

        if self.averaging_on:
            return self.average_viscosity

        return PolymerSolution.rheology_model(shear, self.model["eta_o"], self.model["eta_inf"],
                                      self.model["n"], self.model["l"])  # viscosity [cP]

    def load_rheology_file(self, index: int = 0, is_new: bool=True, override=False):
        # print('PolymerSolution load_rheology_file: index, is_new, override:', index, is_new, override)

        if is_new:
            prev_file = self.file
            super(PolymerSolution, self).set_file()

            if not self.file:
                return -1, False
            elif self.file == prev_file:
                _, _, _, max_index, _ = self.get_rheology_data(0)
                if not override:
                    return max_index, False

        sd, vd, td, max_index, test_name = self.get_rheology_data(index)

        if is_new or override:
            self.initialize_stored(max_index + 1, True)

        gamma = array(sd, dtype='float64')
        self.shear_data = gamma
        self.stored_shear_data[index] = gamma
        viscosity = array(vd, dtype='float64')
        self.viscosity_data = viscosity
        self.stored_viscosity_data[index] = viscosity
        self.ref_temp = td
        self.stored_ref_temps[index] = td
        self.test_name = test_name
        self.stored_test_names[index] = test_name

        try:
            self.fit_viscosity(gamma, viscosity)
        finally:
            pass

        return max_index, is_new

    def get_rheology_data(self, index: int):

        if self.file[-3:] == 'csv':
            return self.get_rheology_data_csv(index)
        else:
            return self.get_rheology_data_excel(index)

    def get_rheology_data_csv(self, index: int):

        file_name = self.file
        rheology = read_csv(file_name, encoding='utf-16', delimiter='\t', header=3)

        s_index = -1
        v_index = -1
        t_index = -1

        for i, key in enumerate(rheology.keys()):
            if key == 'Shear Rate':
                s_index = i
            elif key == 'Viscosity':
                v_index = i
            elif key == 'Temperature':
                t_index = i

        if t_index == -1:
            temp_included = False
        else:
            temp_included = True

        row_index = where(rheology.values[:, 0] == '1')[0]
        if not len(row_index):
            row_index = where(rheology.values[:, 0] == 1.0)[0]

        if index > len(row_index) - 1:
            sd = []
            vd = []
            td = nan
            test = ''
            return sd, vd, td, len(row_index) - 1, test

        current = row_index[index] - 1
        last_row = nan

        if index > 0:
            test = rheology.values[current - 5, 1]
        else:
            header = read_csv(file_name, encoding='utf-16', delimiter='\t', nrows=0)
            test = header.keys()[1]

        for value in rheology.values[row_index[index]:, v_index]:
            current += 1
            try:
                value = float(value)
                print(value)
                if isnan(value):
                    raise ValueError

            except ValueError:
                last_row = current
                break

        if isnan(last_row):
            last_row = current + 1

        sd = rheology.values[row_index[index]:last_row, s_index]
        sd = sd.flatten()
        sd = sd.tolist()
        vd = rheology.values[row_index[index]:last_row, v_index]
        vd = vd.flatten()
        vd = vd.tolist()
        if temp_included:
            td = rheology.values[row_index[index]:last_row, t_index]
            td = td.flatten()
            td = td.tolist()
            for j in range(len(sd)):
                sd[j] = float(sd[j])
                vd[j] = float(vd[j])
                td[j] = float(td[j])
            td = sum(td) / len(td)
        else:
            for j in range(len(sd)):
                sd[j] = float(sd[j])
                vd[j] = float(vd[j])
            td = nan

        return sd, vd, td, len(row_index) - 1, test

    def get_rheology_data_excel(self, index: int):

        rheology = read_excel(self.file)
        v_multiplier = 1.
        is_anton_paar = True
        is_ls1 = False
        if list(rheology.keys())[0] == 'Measuring point':
            # LS-300 data
            is_anton_paar = False
        elif list(rheology.keys())[0] == 'time':
            is_ls1 = True
            is_anton_paar = False

        if is_anton_paar:
            rheology_headers = rheology.values[2, :]
            s_index = where(rheology_headers == 'Shear Rate')
            v_index = where(rheology_headers == 'Viscosity')
            t_index = where(rheology_headers == 'Temperature')
            test = list(rheology.keys())[1]
        elif is_ls1:
            rheology_headers = array(list(rheology.keys()))
            s_index = where(rheology_headers == 'Rate')
            v_index = where(rheology_headers == 'h.1')
            t_index = where(rheology_headers == 'Temp')
            test = ''
        else:
            rheology_headers = array(list(rheology.keys()))
            s_index = where(rheology_headers == 'Shear Rate (1/s)')
            v_index = where(rheology_headers == 'Viscosity (mPas)')
            if v_index[0].size == 0:
                v_index = where(rheology_headers == 'Viscosity (Pas)')
                v_multiplier = 1000.
            t_index = where(rheology_headers == 'Temperature (' + chr(176) + 'C)')
            test = ''

        s_index = s_index[0]
        v_index = v_index[0]
        t_index = t_index[0]

        if len(s_index) == 0:
            return

        if len(v_index) == 0:
            return

        if len(t_index) == 0:
            temp_included = False
        else:
            temp_included = True

        if is_ls1:
            sd = rheology.values[1:, s_index]
            sd = sd.flatten()
            sd = sd.tolist()
            vd = rheology.values[1:, v_index]
            vd = vd.flatten()
            vd = vd.tolist()
            if temp_included:
                td = rheology.values[1:, t_index]
                td = td.flatten()
                td = td.tolist()
                td = sum(td)/len(td)
            else:
                td = nan
            # first_column = rheology.values[:, 0]
            return sd, vd, td, 0, test

        first_column = rheology.values[:, 0]
        row_index = where(first_column == 1.)
        row_index = row_index[0]

        if index > len(row_index) - 1:
            sd = []
            vd = []
            td = nan
            test = ''
            return sd, vd, td, len(row_index) - 1, test

        current = row_index[index] - 1
        last_row = nan
        if index > 0:
            test = rheology.values[current - 5, 1]
        for value in rheology.values[row_index[index]:, 3]:
            current += 1
            if isnan(value) and isnan(last_row):
                last_row = current
                break

        if isnan(last_row):
            last_row = current + 1

        sd = rheology.values[row_index[index]:last_row, s_index]
        sd = sd.flatten()
        sd = sd.tolist()
        vd = rheology.values[row_index[index]:last_row, v_index]
        vd = vd.flatten()
        vd = vd.tolist()
        if temp_included:
            td = rheology.values[row_index[index]:last_row, t_index]
            td = td.flatten()
            td = td.tolist()
            td = sum(td) / len(td)
        else:
            td = nan

        if not is_anton_paar and not is_ls1:
            # For LS-300, need to remove redundant shear rates (both ends) and zero viscosities
            start = 0
            i = 0
            found = False
            while not found:
                if i - 1 < len(sd):
                    if sd[i] != sd[i + 1]:
                        found = True
                        start = i
                    else:
                        i += 1
                else:
                    found = True
                    start = i

            end = len(sd) - 1
            j = len(sd) - 1
            found = False
            while not found:
                if j > 0:
                    if sd[j] != sd[j - 1]:
                        found = True
                        end = j
                    else:
                        j -= 1
                else:
                    found = True
                    end = j

            sd = array(sd[start:end + 1])
            vd = array(vd[start:end + 1])
            no_zeros = where(vd > 0.)
            vd = vd[no_zeros]
            sd = sd[no_zeros]

            sd = list(sd)
            vd = v_multiplier * vd
            vd = list(vd)

        return sd, vd, td, len(row_index) - 1, test

    def fit_viscosity(self, gamma, viscosity):

        if not self.averaging_on:
            fit_params, _ = curve_fit(PolymerSolution.rheology_model, gamma, viscosity,
                                      bounds=(0, [10000., 100., 1., 100.]))
            self.model["eta_o"] = fit_params[0]
            self.model["eta_inf"] = fit_params[1]
            self.model["n"] = fit_params[2]
            self.model["l"] = fit_params[3]
        else:
            self.model["eta_o"] = nan
            self.model["eta_inf"] = nan
            self.model["n"] = nan
            self.model["l"] = nan

        self.average_viscosity = average(viscosity)

    def plot_viscosity(self):
        gamma = logspace(-1, 3)
        fit_viscosity = PolymerSolution.rheology_model(gamma, self.model["eta_o"], self.model["eta_inf"],
                                               self.model["n"], self.model["l"])
        plt.plot(gamma, fit_viscosity)
        # plt.plot(gamma, viscosity)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Shear Rate [1/s]')
        plt.ylabel('Viscosity [cP]')
        plt.show()

    @staticmethod
    def rheology_model(shear, eta_o, eta_inf, n, l):
        return eta_inf + divide((eta_o - eta_inf), (1. + power(l * shear, n)))

    @staticmethod
    def temp_dependence_polymer(temp, temp_o):
        return divide(temp_o, temp)

    @staticmethod
    def temp_dependence_oil(temp, m, b):
        return power(10., power(10., m * temp + b))


class PolymerSolutionView(QMainWindow, utils.SignalView):
    """This class creates the main GUI as a subclass of QMainWindow. The central widget contains QAxes
    for displaying the data, back and next QButtons for traversing multiple measurements in a data file,
    a submit button for submitting information to the database, and a ShearViscosityWidget for entering
    data."""

    # Note: should consider augmenting PolymerView to include the interactive QChart, prev and next
    # QButtons, and likely the polymer and brine name text, probably as edit boxes.

    # Note: this needs temperature edit for files without that data (or manually entered). It also
    # needs a concentration edit. Both of these could go in PolymerView (see note above).

    # Note: a QGridLayout would help for the central widget.

    # Note: database_action should probably be moved out of the class and into the main function in the
    # module. The operation would be more transparent.

    # Note: the ShearViscosityTable should be able to spit out its own data as is queried in
    # update_data_from_table.

    # Note: in update_data_from_plot, why discriminate against an empty table?

    def __init__(self, rheo_fluid, parent=None):
        """The constructor creates the UI without making it visible so that the user can be asked
        whether they want to start from an existing database or create a new one."""

        super(PolymerSolutionView, self).__init__(parent=parent)

        #################################
        ########## BEGIN SETUP ##########
        #################################

        # Basic window properties (size, title, icon)

        # print('PolymerSolutionView parent:', parent)
        try:
            if isinstance(parent, PolymerRheologyView):
                sf = parent.width() / 200.
            elif isinstance(parent.parent(), OilSampleTool):
                sf = parent.parent().sf
            else:
                sf = parent.parent().width() / 280.
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()
            try:
                sf = parent.width() / 900.
            except Exception as e2:
                QMessageBox(parent=self, text=str(e2)).exec_()
                sf = 1.

        self.setFixedSize(int(sf * 550), int(sf * 550))
        self.setWindowIcon(QIcon('Coreholder Cropped.jpg'))

        ##### Menus #####

        # file menu with new and previous database actions
        # self.file_menu = QMenu('&File')
        # self.file_menu.addAction('&New DB', self.new_database)
        # self.file_menu.addAction('&Prev DB', self.prev_database)
        # self.menuBar().addMenu(self.file_menu)

        # tools menu with rheology import, brine entry/selection, and polymer entry/selection actions
        # self.tools_menu = QMenu('&Tools')
        # self.tools_menu.addAction('&Clear UI...', self.clear_ui)
        # self.tools_menu.addAction('&Import Rheo...', self.load_rheology_file)
        # self.tools_menu.addAction('&Brines...', self.brine_entry)
        # self.tools_menu.addAction('&Polymers...', self.polymer_entry)
        # self.menuBar().addMenu(self.tools_menu)

        ##### The Central Widget #####

        # chart, view, and axes with titles
        self.chart = QChart(flags=Qt.WindowFlags())     # chart
        self.view = QChartView(self.chart)              # view
        self.x_axis = QLogValueAxis()                   # axes
        self.y_axis = QLogValueAxis()
        self.chart.setAxisX(self.x_axis)
        self.chart.setAxisY(self.y_axis)
        self.chart.legend().hide()                      # no legend
        self.x_axis.setTitleText(chr(947) + u' [1/s]')  # axes titles
        self.y_axis.setTitleText(chr(951) + u' [cP]')
        f = self.y_axis.titleFont()                     # title font size
        f.setPointSize(12)
        self.x_axis.setTitleFont(f)
        self.y_axis.setTitleFont(f)

        # assembling the the widget with chart, brine and polymer info, and submit button
        left_widget = QWidget()                 # widget to hold everything on left side of GUI
        left_layout = QVBoxLayout()             # layout for left widget is vertical
        left_widget.setLayout(left_layout)      # add set left layout as layout for left widget

        label_box = QGroupBox(parent=self)      # group box for brine and polymer info
        if isinstance(rheo_fluid.base_fluid, Brine):
            self.setWindowTitle('Polymer Solution Rheology Tool')
            self.base_fluid_label = QLabel(parent=self, text='Brine: ' + rheo_fluid.base_fluid.name)
            if rheo_fluid.primary_additive is not None:
                self.additive_label = QLabel(parent=self, text='Polymer: ' + rheo_fluid.primary_additive.name)
            else:
                self.additive_label = QLabel(parent=self, text='Polymer: None')
        else:
            self.setWindowTitle('Oil Rheology Tool')
            self.base_fluid_label = QLabel(parent=self, text='Oil: ' + rheo_fluid.base_fluid.name)
            if rheo_fluid.primary_additive is not None:
                self.additive_label = QLabel(parent=self, text='Diluent: ' + rheo_fluid.primary_additive.name)
            else:
                self.additive_label = QLabel(parent=self, text='Diluent: None')
        label_layout = QHBoxLayout()            # horizontal layout for brine and polymer info
        label_layout.addWidget(self.base_fluid_label)    # add brine and polymer info to label layout
        label_layout.addWidget(self.additive_label)
        label_box.setLayout(label_layout)           # add label layout to label group box
        left_layout.addWidget(label_box)            # add label group box to left widget
        left_layout.addWidget(self.view)            # add chart view to left widget layout
        self.submit_button = QPushButton(parent=self, text='Import Rheology...')    # submit button
        self.submit_button.clicked.connect(lambda: self.load_rheology_file(True))     # submit button slot connection
        self.tool_parent = None
        self.oil_if = None
        left_layout.addWidget(self.submit_button)   # add submit button to left widget layout

        layout = QHBoxLayout()  # main (overall) layout is horizontal

        self.prev_button = QPushButton(parent=self, text='<')   # previous data button w/ slot connect
        self.prev_button.setFixedWidth(15)
        self.prev_button.clicked.connect(self.load_rheology_file_wrapper_prev)
        layout.addWidget(self.prev_button)

        layout.addWidget(left_widget)   # add the left widget to the layout

        self.next_button = QPushButton(parent=self, text='>')   # next data button w/ slot connect
        self.next_button.setFixedWidth(15)
        self.next_button.clicked.connect(self.load_rheology_file_wrapper_next)
        layout.addWidget(self.next_button)

        # self.table_widget = ShearViscosityWidget(self)      # the ShearViscosityWidget for shear data
        # self.table_widget.setFixedSize(QSize(250, 500))
        # layout.addWidget(self.table_widget)

        c_widget = QWidget()                # at last, the main widget associated with the main layout
        c_widget.setLayout(layout)
        self.setCentralWidget(c_widget)     # set as central widget for QMainWindow

        # The temperature and concentration labels/edits are super-imposed over the top-right corner of
        # the QAxes. They are connected to slots that modify polymer ref_temp and concentration.
        self.temp_label = QLabel(parent=self, text=u'Temp [' + chr(176) + 'C]:')
        self.temp_edit = QLineEdit(parent=self)
        self.temp_edit.textChanged.connect(self.update_temperature)
        self.temp_label.setFixedSize(int(sf * 60), int(sf * 20))
        self.temp_edit.setFixedSize(int(sf * 30), int(sf * 20))
        self.temp_label.move(int(sf * 400), int(sf * 90))
        self.temp_edit.move(int(sf * 465), int(sf * 90))

        # self.conc_label = QLabel(parent=self, text=u'Conc. [ppm]:')
        # self.conc_edit = QLineEdit(parent=self)
        # self.conc_edit.textChanged.connect(self.update_concentration)
        # self.conc_label.setFixedSize(60, 20)
        # self.conc_edit.setFixedSize(30, 20)
        # self.conc_label.move(400, 115)
        # self.conc_edit.move(465, 115)

        self.averaging_checkbox = QCheckBox(parent=self, text='Use Average?')     # designate sheared polymer
        self.averaging_checkbox.move(int(sf * 430), int(sf * 40))
        self.averaging_checkbox.setChecked(rheo_fluid.averaging_on)
        self.averaging_checkbox.clicked.connect(lambda: self.toggle_averaging(False))

        # Create a display for the model fit parameters.
        fit_params_string = chr(951) + u'_o: cP ;  ' + \
                            chr(951) + u'_inf: cP ;  ' + u'n:  ;  l: '
        self.fit_params_label = QLabel(parent=self, text=fit_params_string)
        self.fit_params_label.setFixedWidth(int(sf * 300))
        self.fit_params_label.move(int(sf * 50), int(sf * 80))

        # Create display text to let the user know whether the data is already in the database.
        # self.already_entered_label = QLabel(parent=self, text='')
        # self.already_entered_label.setStyleSheet('background: white')
        # self.already_entered_label.setAlignment(Qt.AlignCenter)
        # self.already_entered_label.setFixedSize(135, 20)
        # self.already_entered_label.move(55, 450)

        ##### Other objects and quantities to track #####

        self.database_name = ''  # the name of the database file, to be displayed on window title
        self.directory = ''  # the directory of the database
        self.database = None  # the handle to the database connection
        self.allow_data_exclusion = True

        self.polymer_solution = rheo_fluid  # the Polymer object that imports, stores, and fits the data
        # self.polymer_solution.fit_changed.connect(self.update_fit_params_label)
        # self.update_fit_params_label()

        self.view_mode = False
        if self.polymer_solution.stored_shear_data:
            self.view_mode = True
            self.n_entries = len(self.polymer_solution.stored_shear_data) - 1
            self.current_entry = 0
            self.polymer_solution.average_viscosity = self.polymer_solution.stored_averages[0]
            self.polymer_solution.model = copy(self.polymer_solution.stored_models[0])
            self.polymer_solution.ref_temp = self.polymer_solution.stored_ref_temps[0]
            self.polymer_solution.shear_data = self.polymer_solution.stored_shear_data[0]
            self.polymer_solution.viscosity_data = self.polymer_solution.stored_viscosity_data[0]
            self.polymer_solution.test_name = self.polymer_solution.stored_test_names[0]
            self.polymer_solution.excluded = self.polymer_solution.stored_excluded[0]
        else:
            self.n_entries = -1  # the number of measurements in the AntonPaar file
            self.current_entry = -1  # the measurement number in the file being displayed
        # self.brine_name = ''  # the name of the brine, to be displayed above the QAxes
        # self.polymer_name = ''  # the name of the polymer, to be displayed above the QAxes
        # self.brine_id = 0  # the id of the brine, referenced in the database 'brines' table
        # self.polymer_id = 0  # the id of the polymer, referenced in the database 'polymers' table
        self.plot_curves = []  # the list that will hold handles to the datapoint QScatterSeries
        self.fit_curve = None  # the placeholder for the fit line QLineSeries

        # self.brine_tool = None  # the placeholder for the DatabaseBrineTool
        # self.polymer_tool = None  # the placeholder for the DatabasePolymerTool

        try:
            print(self.polymer_solution.ref_temp)
        except Exception as e:
            QMessageBox(parent=self, text='No ref temp: ' + str(e)).exec_()

        try:
            print(self.polymer_solution.file)
        except Exception as e:
            QMessageBox(parent=self, text='No file: ' + str(e)).exec_()

        try:
            print(self.polymer_solution.shear_data)
        except Exception as e:
            QMessageBox(parent=self, text='No shear data: ' + str(e)).exec_()

        if not isnan(self.polymer_solution.ref_temp):
            self.temp_edit.setText(str(np_round(self.polymer_solution.ref_temp, 1)))  # have temp and conc edits reflect polymer
        # self.conc_edit.setText(str(int(self.polymer_solution.concentration)))
        # self.conc_edit.setReadOnly(True)

        if self.polymer_solution.file:
            self.show_file_name()   # show the empty string

        if isinstance(self.polymer_solution.shear_data, ndarray):
            self.set_plot_data()

        #################################
        ########### END SETUP ###########
        #################################

        # show the polymer rheology tool
        self.show()
        # print('made it to end of PolymerSolutionView setup.')

    def update_temperature(self):
        """This function updates when the text in the temp_edit is changed. It tries to convert the
        text string into an float and then checks that the float is non-negative before permitting it.
        Otherwise, the value is reset to what it was prior."""

        is_valid = False    # initialize valid entry boolean marker

        try:
            temperature = float(self.temp_edit.text())  # try interpreting temp_edit text as float
            if temperature >= 0.:
                # If the temperature is non-negative in Celsius, record as the polymer ref_temp
                self.polymer_solution.ref_temp = temperature
                is_valid = True
            # print(self.polymer.ref_temp)    # print the new polymer ref_temp
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()
            # print(self.temp_edit.text(), e)     # report any error if Exception

        # Reset the temp_edit text if the entry was invalid.
        if not is_valid:
            self.temp_edit.blockSignals(True)   # block signals for special case of NaN initialization
            self.temp_edit.setText(str(self.polymer.ref_temp))
            self.temp_edit.blockSignals(False)

    def show_file_name(self):
        """This function simply shows the user the name of the file being imported."""
        i = self.polymer_solution.file.rfind('/')   # polymer.file is a full file path -- find last '/' before file name
        # send file and test name as message
        self.statusBar().showMessage(self.polymer_solution.file[i + 1:] + ':   ' +
                                     self.polymer_solution.stored_test_names[self.current_entry])

    def update_fit_params_label(self):
        """This function is called when a signal is generated from the Polymer object (self.polymer)
        indicating that fit_changed."""

        m = self.polymer_solution.model  # the dictionary of model parameters
        if not isnan(m['eta_o']):
            v = (m['eta_o'], m['eta_inf'], m['n'], m['l'])  # construct tuple for string formatting
            fit_params_string = chr(951) + '_o: {0:7.2f} cP ;  ' + chr(951) + \
                '_inf: {1:7.2f} cP ;  n: {2:7.4f} ;  l: {3:7.4f}'   # string with formatting
        else:
            v = (self.polymer_solution.average_viscosity, )
            fit_params_string = chr(951) + ': {0:7.2f} cP'

        self.fit_params_label.setText(fit_params_string.format(*v))  # display formatted string

    def toggle_averaging(self, silent: bool=True):

        # print('toggling averaging.')

        if self.averaging_checkbox.isChecked():
            # print('setting averaging to ON.')
            self.polymer_solution.averaging_on = True
        else:
            # print('setting averaging to OFF.')
            self.polymer_solution.averaging_on = False

        if not silent:
            self.set_plot_data()

    def load_rheology_file_wrapper_prev(self):
        """This is the slot for the previous data entry button. It is a wrapper for load_rheology_file."""

        if self.n_entries == -1:
            # If there's no data, leave the function.
            return

        if self.current_entry < 1:
            # If the current entry is the first, don't try to decrement.
            return

        if not isnan(self.polymer_solution.stored_averages[self.current_entry - 1]) \
                and isnan(self.polymer_solution.stored_models[self.current_entry - 1]['eta_o']):
            self.averaging_checkbox.setChecked(True)
        elif not isnan(self.polymer_solution.stored_models[self.current_entry - 1]['eta_o']):
            self.averaging_checkbox.setChecked(False)
        self.toggle_averaging()

        # Call load_rheology_file with is_new=False and the entry as one less than the current.
        self.load_rheology_file(False, self.current_entry - 1)

    def load_rheology_file_wrapper_next(self):
        """This is the slot for the next data entry button. It is a wrapper for load_rheology_file."""

        if self.n_entries == -1:
            # If there's no data, leave the function.
            return

        if self.current_entry >= self.n_entries:
            # If the current entry is the last, don't try to increment.
            return

        elif isnan(self.polymer_solution.stored_averages[self.current_entry + 1]) and \
                isnan(self.polymer_solution.stored_models[self.current_entry + 1]['eta_o']) and self.view_mode:
            return

        # print(self.polymer_solution.stored_averages, self.current_entry, self.polymer_solution.stored_models)
        if not isnan(self.polymer_solution.stored_averages[self.current_entry + 1]) \
                and isnan(self.polymer_solution.stored_models[self.current_entry + 1]['eta_o']):
            self.averaging_checkbox.setChecked(True)
        elif not isnan(self.polymer_solution.stored_models[self.current_entry + 1]['eta_o']):
            self.averaging_checkbox.setChecked(False)

        self.toggle_averaging()

        # Call load_rheology_file with is_new=False and the entry as one more than the current.
        self.load_rheology_file(False, self.current_entry + 1)

    def load_rheology_file(self, is_new: bool=True, entry: int=0):
        """This function calls the load_rheology_file method of the polymer object and sets the data
        displayed on the plot and table."""

        # This makes sure the data is not out-of-range.
        # Note: may be unnecessary. Is this called somewhere else with is_new=False?
        view_mode_toggled = False
        if not is_new:
            if entry > self.n_entries:
                return
        else:
            if self.view_mode:
                self.view_mode = False
                view_mode_toggled = True

        new_loaded = 0
        n_entries = 0
        if not self.view_mode:
            try:
                # actually load the file and save n_entries
                n_entries, new_loaded = \
                    self.polymer_solution.load_rheology_file(entry, is_new, override=view_mode_toggled)

            except Exception as e:
                QMessageBox(parent=self, text=str(e)).exec_()
                return

        if is_new:
            if not new_loaded and not view_mode_toggled:
                # If the user canceled out of the file dialog, there's nothing to do.
                return
            else:
                # If the file is new (and was loaded), then n_entries (the number of entries in the
                # file) needs to be updated.
                self.n_entries = n_entries

        self.current_entry = entry  # save the current entry

        self.show_file_name()  # show the name of the new file

        self.polymer_solution.excluded = self.polymer_solution.stored_excluded[entry]

        self.set_plot_data()        # set the plot data to the new data

        # Set the temp_edit to reflect polymer ref_temp.
        self.temp_edit.blockSignals(True)
        self.temp_edit.setText(str(np_round(self.polymer_solution.stored_ref_temps[self.current_entry], 1)))
        self.temp_edit.blockSignals(False)

    def rheo_selection_tool_launcher(self):

        OilInjectionFluidRheologySelectionTool(self.tool_parent, self.oil_if, self)

    def copy_rheology(self, rheo_fluid: PolymerSolution):

        self.polymer_solution.average_viscosity = rheo_fluid.average_viscosity
        self.polymer_solution.stored_averages = deepcopy(rheo_fluid.stored_averages)
        self.polymer_solution.concentration = rheo_fluid.concentration
        self.polymer_solution.model = deepcopy(rheo_fluid.model)
        self.polymer_solution.stored_models = deepcopy(rheo_fluid.stored_models)
        self.polymer_solution.ref_temp = rheo_fluid.ref_temp
        self.polymer_solution.stored_ref_temps = deepcopy(rheo_fluid.stored_ref_temps)
        self.polymer_solution.shear_data = deepcopy(rheo_fluid.shear_data)
        self.polymer_solution.stored_shear_data = deepcopy(rheo_fluid.stored_shear_data)
        self.polymer_solution.viscosity_data = deepcopy(rheo_fluid.viscosity_data)
        self.polymer_solution.stored_viscosity_data = deepcopy(rheo_fluid.stored_viscosity_data)
        self.polymer_solution.test_name = copy(rheo_fluid.test_name)
        self.polymer_solution.stored_test_names = copy(rheo_fluid.stored_test_names)
        self.polymer_solution.excluded = copy(rheo_fluid.excluded)
        self.polymer_solution.stored_excluded = copy(rheo_fluid.stored_excluded)

        ps_view = PolymerSolutionView(rheo_fluid=self.polymer_solution, parent=self.parent())
        ps_view.submit_button.clicked.disconnect()
        ps_view.tool_parent = self.tool_parent
        ps_view.oil_if = self.oil_if
        ps_view.submit_button.clicked.connect(ps_view.rheo_selection_tool_launcher)
        ps_view.allow_data_exclusion = False
        self.close()

    def set_plot_data(self):
        """This sets the plot data based on the data in the Polymer object. It does some axis ranging."""

        # Remove all series from the chart and clear references to those series.
        self.chart.removeAllSeries()
        self.plot_curves = []
        self.fit_curve = None

        # pulling the polymer shear and viscosity data
        vd = self.polymer_solution.stored_viscosity_data[self.current_entry]
        sd = self.polymer_solution.stored_shear_data[self.current_entry]

        if isinstance(vd, float) or isinstance(sd, float):
            # return if either viscosity or shear data is NaN
            return

        if sd.size != vd.size:
            # return if viscosity or shear data aren't the same size
            return

        k = 0
        for x, y in zip(sd, vd):
            # Loop through the data one point at a time, and add single-point scatter curves to the plot.
            curve = QScatterSeries()
            curve.setUseOpenGL(False)   # OpenGL fails for setting marker size
            if k in self.polymer_solution.excluded:
                curve.setColor(Qt.lightGray)
            else:
                curve.setColor(Qt.black)    # the data are black unless they are graphically excluded
            curve.setMarkerSize(10.)    # big markers
            curve.setMarkerShape(QScatterSeries.MarkerShapeCircle)  # the markers are circles
            self.chart.addSeries(curve)     # add the single-point series to the axes
            curve.append(QPointF(x, y))
            curve.attachAxis(self.x_axis)   # attach the axes
            curve.attachAxis(self.y_axis)
            self.plot_curves.append(curve)  # keep a reference to the data curve for user interactions
            k += 1

        # Use log2 axis limit definitions for this log-log plot.
        base = 2.

        # Each limit is one order of two above or below the max or min of the data.
        order_y_max = utils.find_order_above(vd.max(), base=base)
        self.y_axis.setMax(base ** order_y_max)

        order_y_min = utils.find_order_below(vd.min(), base=base)
        self.y_axis.setMin(base ** order_y_min)

        order_x_max = utils.find_order_above(sd.max(), base=base)
        self.x_axis.setMax(base ** order_x_max)

        order_x_min = utils.find_order_below(sd.min(), base=base)
        self.x_axis.setMin(base ** order_x_min)

        if sd.size < 4 and self.polymer_solution.averaging_on:
            # If there are fewer than four data points, don't try a fit with Carreau model.
            pass
        else:
            # Fit the data and add the fit curve.
            try:
                exclude = self.polymer_solution.excluded
                # remaining points to be fitted
                fit_x = delete(self.polymer_solution.stored_shear_data[self.current_entry], exclude)
                fit_y = delete(self.polymer_solution.stored_viscosity_data[self.current_entry], exclude)
                self.polymer_solution.fit_viscosity(fit_x, fit_y)
                self.polymer_solution.stored_excluded[self.current_entry] = self.polymer_solution.excluded
                self.polymer_solution.stored_averages[self.current_entry] = self.polymer_solution.average_viscosity

                # if not self.polymer_solution.averaging_on:
                self.polymer_solution.stored_models[self.current_entry] = copy(self.polymer_solution.model)

                self.plot_fit_curve()
            except Exception as e:
                QMessageBox(parent=self, text=str(e)).exec_()
                pass

    def plot_fit_curve(self, base: float=2.):
        """This function generates a high-resolution curve that plots the model fit to the polymer
        viscosity data."""

        if self.fit_curve is not None:
            # Start by removing the previous fit curve, if there is one.
            self.fit_curve.detachAxis(self.x_axis)
            self.fit_curve.detachAxis(self.y_axis)
            self.chart.removeSeries(self.fit_curve)

        # Find bounds for a log space based on current axis limits, defined in set_plot_data by powers
        # of 2 (hence the default value of 'base').
        order_x_min = log(self.x_axis.min(), base)
        order_x_max = log(self.x_axis.max(), base)
        x = logspace(order_x_min, order_x_max, base=base)

        # Pull the model parameters from the polymer rheology model, and generate viscosities pertaining
        # to the shear rates in 'x'.
        params = self.polymer_solution.model
        if not isnan(params['eta_o']):
            y = self.polymer_solution.rheology_model(x, params["eta_o"], params["eta_inf"], params["n"], params["l"])
        else:
            y = self.polymer_solution.stored_averages[self.current_entry] * ones(x.shape)
        fit_polyline = utils.series_to_polyline(x, y)   # create a polyline to load into the curve

        # Define a fit curve as a QLineSeries
        fit_curve = QLineSeries()
        pen = fit_curve.pen()
        pen.setColor(Qt.red)    # the curve should be red
        pen.setWidthF(1.)       # and rather thin
        fit_curve.setPen(pen)
        fit_curve.setUseOpenGL(False)   # OpenGL causes problems sometimes
        fit_curve.append(fit_polyline)  # load the modeled data

        self.chart.addSeries(fit_curve)     # add the fit curve to the chart and attach the axes
        fit_curve.attachAxis(self.x_axis)
        fit_curve.attachAxis(self.y_axis)

        self.fit_curve = fit_curve  # store a reference to the fit curve
        self.update_fit_params_label()

    def mousePressEvent(self, a0: QtGui.QMouseEvent):
        """This function handles manual exclusion of shear viscosity data by the user, initiated by
        clicking near points on the plot. The color of individual points is toggled from black to gray
        and back as points are excluded or included."""

        if not self.polymer_solution.file:
            return

        if not self.allow_data_exclusion:
            return

        self.show_file_name()   # this is an easy way to re-print the file name if it is lost

        x_vector = log10(self.polymer_solution.shear_data)
        y_vector = log10(self.polymer_solution.viscosity_data)

        if type(x_vector) == np_float:
            return  # return if no data

        x = a0.x() - 40     # x and y points where user clicked, with offset adjustment
        y = a0.y() - 80
        point = QPointF(x, y)
        point = self.view.mapFromScene(point)   # map the mouse point to graph coordinates
        point = self.chart.mapToScene(point)
        point = self.chart.mapToValue(point)

        if point.x() < self.x_axis.min() or point.x() > self.x_axis.max():
            return  # return if outside x-axis bounds

        if point.y() < self.y_axis.min() or point.y() > self.y_axis.max():
            return  # return if outside y-axis bounds

        x = log10(point.x())     # transform x and y to log10 scale for distance computation
        y = log10(point.y())

        index = utils.find_nearest(x_vector, y_vector, x, y)    # find the point that will be toggled
        data_point = self.plot_curves[index]
        if data_point.color() == Qt.black:
            data_point.setColor(Qt.lightGray)
        else:
            data_point.setColor(Qt.black)

        exclude = []    # generate a list of the excluded points for the fit based on gray color
        for i, curve in enumerate(self.plot_curves):
            if curve.color() == Qt.lightGray:
                exclude.append(i)

        fit_x = delete(self.polymer_solution.shear_data, exclude)     # remaining points to be fitted
        fit_y = delete(self.polymer_solution.viscosity_data, exclude)

        if (fit_x.size < 4 and not self.polymer_solution.averaging_on) or fit_x.size == 0:
            # Don't fit if fewer than four points remaining.
            msg = QMessageBox(parent=self, text='Too few points.')
            msg.show()
            # Reset the data point.
            if data_point.color() == Qt.black:
                data_point.setColor(Qt.lightGray)
            else:
                data_point.setColor(Qt.black)

            return  # return without fitting if fewer than 4 points

        try:
            self.polymer_solution.fit_viscosity(fit_x, fit_y)    # try fitting the remaining points
            self.polymer_solution.stored_averages[self.current_entry] = self.polymer_solution.average_viscosity
            if not self.polymer_solution.averaging_on:
                self.polymer_solution.stored_models[self.current_entry] = copy(self.polymer_solution.model)
        except Exception as e:
            # If the fit fails, report the error to the user as a message.
            err_string = 'Fitting failed!<br><br>' + str(e)
            msg = QMessageBox(parent=self, text=err_string)
            msg.setWindowTitle('Error')
            msg.show()
            # Reset the data point the user clicked on.
            if data_point.color() == Qt.black:
                data_point.setColor(Qt.lightGray)
            else:
                data_point.setColor(Qt.black)
            return

        # Plot the fit curve. If "try" is successful, this will plot the new curve. If "except"
        # executes, the same curve will be re-plotted.
        self.plot_fit_curve()
        self.polymer_solution.excluded = exclude
        self.polymer_solution.stored_excluded[self.current_entry] = exclude

    def submit(self):
        """This function executes when the submit button is pressed. It composes a SQLite INSERT query
        that sends a full row entry to the viscosity_measurements table."""

        print(self)

    def closeEvent(self, a0: QtGui.QCloseEvent):
        print(self, a0)
        self.close_signal.emit(self)


class OilInjectionFluidRheologySelectionTool(QDialog):

    def __init__(self, parent: OilInjectionFluidView, oil_if: InjectionFluid, ps_view: PolymerSolutionView):
        super(OilInjectionFluidRheologySelectionTool, self).__init__(parent=parent)
        self.ps_view = ps_view
        self.oil_if = oil_if

        self.setLayout(QVBoxLayout())
        self.setWindowTitle('Select Rheo...')

        txt = QLabel(parent=self, text=oil_if.specific_fluid.oil_sample.name + ' Rheology List:')
        self.listbox = QListWidget(parent=self)
        rheologies_names = oil_if.specific_fluid.oil_sample.rheologies_list.signal_lists[0].item_names
        self.listbox.addItems(rheologies_names)
        self.listbox.itemDoubleClicked.connect(self.rheology_selected)

        self.layout().addWidget(txt)
        self.layout().addWidget(self.listbox)

        self.setFixedSize(int(250 * parent.sf), int(250 * parent.sf))

        self.show()

    def rheology_selected(self, item: QListWidgetItem):

        rheologies_names = self.oil_if.specific_fluid.oil_sample.rheologies_list.signal_lists[0].item_names
        rheologies_list = self.oil_if.specific_fluid.oil_sample.rheologies_list.signal_lists[0].objects
        i = rheologies_names.index(item.text())
        rheology = rheologies_list[i]

        self.ps_view.copy_rheology(rheo_fluid=rheology.specific_fluid.polymer_solution)
        self.close()


class PolymerRheologyImportTool(PolymerSolutionView):

    def __init__(self, rheo_fluid, parent=None, polymer_rheology=None):
        self.stored_concentrations = []
        # print('PolymerRheologyImportTool parent:', parent)
        if polymer_rheology is not None:
            if polymer_rheology.polymer_solutions_list:
                rheo_fluid = PolymerRheologyImportTool.create_single_solution(polymer_rheology)
                for solution in polymer_rheology.polymer_solutions_list:
                    for _ in solution.stored_models:
                        self.stored_concentrations.append(solution.concentration)
        self.stored_concentrations = array(self.stored_concentrations)
        super(PolymerRheologyImportTool, self).__init__(rheo_fluid=rheo_fluid, parent=parent)
        self.polymer_rheology = polymer_rheology

        self.setWindowTitle('Polymer Rheology View')
        sf = parent.width() / 200.
        self.setFixedSize(int(sf * 550), int(sf * 550))
        self.conc_label = QLabel(parent=self, text=u'Conc. [ppm]:')
        self.conc_edit = QLineEdit(parent=self)
        self.conc_label.setFixedSize(int(sf * 60), int(sf * 20))
        self.conc_edit.setFixedSize(int(sf * 30), int(sf * 20))
        self.conc_edit.editingFinished.connect(self.update_concentration)
        self.layout().addWidget(self.conc_label)
        self.layout().addWidget(self.conc_edit)
        self.conc_label.move(int(sf * 400), int(sf * 115))
        self.conc_edit.move(int(sf * 465), int(sf * 115))

        self.update_concentration_from_stored()

    def load_rheology_file(self, is_new: bool=True, entry: int=0):

        ce_before = self.current_entry
        super(PolymerRheologyImportTool, self).load_rheology_file(is_new, entry)

        if self.current_entry == ce_before and not is_new:
            return

        if is_new:
            if self.n_entries > -1:
                self.stored_concentrations = zeros(self.n_entries + 1)
                self.stored_concentrations[:] = nan

        try:
            self.update_concentration_from_stored()
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()

    def update_concentration(self):

        txt = self.conc_edit.text()

        if not txt:
            self.stored_concentrations[self.current_entry] = nan
            return

        try:
            val = float(txt)
            val = int(val)
            if val < 0:
                raise ValueError

            print('stored concentrations:', self.stored_concentrations, 'val:', val)
            self.stored_concentrations[self.current_entry] = val

        except ValueError as e:
            print(e)
            self.update_concentration_from_stored()

    def update_concentration_from_stored(self):

        if not self.stored_concentrations.tolist():
            return

        stored = self.stored_concentrations[self.current_entry]
        self.conc_edit.blockSignals(True)
        if isnan(stored):
            self.conc_edit.setText('')
        else:
            self.conc_edit.setText(str(int(stored)))
        self.conc_edit.blockSignals(False)

    @staticmethod
    def create_solutions_list(stored_concentrations: array, polymer_solution: PolymerSolution, poly: Polymer):
        non_nans = where(~isnan(stored_concentrations))[0]

        if not non_nans.size:
            # self.polymer_rheology.polymer_solutions_list = []
            # print(self.polymer_rheology.polymer_solutions_list)
            return []

        scs = stored_concentrations[non_nans]
        sts = array(polymer_solution.stored_ref_temps)
        sts = sts[non_nans]
        savs = array(polymer_solution.stored_averages)
        savs = savs[non_nans]
        sms = []
        ssd = []
        svd = []
        stns = []
        se = []

        for non_nan in non_nans:
            sms.append(polymer_solution.stored_models[non_nan])
            ssd.append(polymer_solution.stored_shear_data[non_nan])
            svd.append(polymer_solution.stored_viscosity_data[non_nan])
            stns.append(polymer_solution.stored_test_names[non_nan])
            se.append(polymer_solution.stored_excluded[non_nan])

        uscs = unique(scs)
        psl = []

        for usc in uscs:
            psl.append(PolymerSolution(base_fluid=polymer_solution.base_fluid,
                                       primary_additive=poly, concentration=usc))

            same_c = where(scs == usc)[0]
            new_models = []
            new_shear_data = []
            new_viscosity_data = []
            new_test_names = []
            new_stored_excluded = []

            for ind in same_c:
                new_models.append(sms[ind])
                new_shear_data.append(ssd[ind])
                new_viscosity_data.append(svd[ind])
                new_test_names.append(stns[ind])
                new_stored_excluded.append(se[ind])

            psl[-1].stored_models = new_models
            psl[-1].stored_averages = savs[same_c]
            psl[-1].stored_ref_temps = sts[same_c]
            psl[-1].stored_shear_data = new_shear_data
            psl[-1].stored_viscosity_data = new_viscosity_data
            psl[-1].stored_test_names = new_test_names
            psl[-1].stored_excluded = new_stored_excluded

            psl[-1].name = polymer_solution.name
            psl[-1].dir = polymer_solution.dir
            psl[-1].file = polymer_solution.file

        # for ps in psl:
        #     print(ps.concentration, ps.stored_ref_temps, ps.stored_models, ps.stored_averages, ps.stored_shear_data,
        #           ps.stored_viscosity_data, ps.stored_test_names, ps.stored_excluded)

        return psl

    @staticmethod
    def create_single_solution(polymer_rheology: PolymerRheology):
        ps = PolymerSolution(base_fluid=polymer_rheology.base_fluid.specific_fluid.brine,
                             primary_additive=polymer_rheology.poly)

        sts = []
        savs = []
        sms = []
        ssd = []
        svd = []
        stns = []
        se = []

        for solution in polymer_rheology.polymer_solutions_list:
            for ref_temp, average, model, shear_datum, viscosity_datum, test_name, excluded in \
                    zip(solution.stored_ref_temps, solution.stored_averages, solution.stored_models,
                        solution.stored_shear_data, solution.stored_viscosity_data, solution.stored_test_names,
                        solution.stored_excluded):
                sts.append(ref_temp)
                savs.append(average)
                sms.append(model)
                ssd.append(array(shear_datum))
                svd.append(array(viscosity_datum))
                stns.append(test_name)
                se.append(excluded)

        ps.stored_ref_temps = sts
        ps.stored_averages = savs
        ps.stored_models = sms
        ps.stored_shear_data = ssd
        ps.stored_viscosity_data = svd
        ps.stored_test_names = stns
        ps.stored_excluded = se

        ps.name = polymer_rheology.polymer_solutions_list[0].name
        ps.dir = polymer_rheology.polymer_solutions_list[0].dir
        ps.file = polymer_rheology.polymer_solutions_list[0].file

        # print(ps.stored_ref_temps, ps.stored_averages, ps.stored_models, ps.stored_shear_data,
        #       ps.stored_viscosity_data, ps.stored_test_names, ps.stored_excluded)

        return ps

    def closeEvent(self, a0: QtGui.QCloseEvent):
        super(PolymerRheologyImportTool, self).closeEvent(a0)

        self.polymer_rheology.polymer_solutions_list = \
            PolymerRheologyImportTool.create_solutions_list(self.stored_concentrations, self.polymer_solution,
                                                            self.polymer_rheology.poly)

        try:
            self.parent().refresh_viscosity_estimation()
            self.parent().import_tool = None
            if len(self.polymer_rheology.polymer_solutions_list) > 1:
                self.parent().view_button.setEnabled(True)
            else:
                self.parent().view_button.setEnabled(False)
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()


if __name__ == '__main__':
    app = QApplication(argv)
    b = Brine('DI')
    polymer = Polymer('FP3630S')
    ps = PolymerSolution(brine=b, polymer=polymer, concentration=1000)
    psv = PolymerSolutionView(ps)
    sys_exit(app.exec_())
