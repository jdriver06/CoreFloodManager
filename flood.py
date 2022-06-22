
from __future__ import division
from PyQt5.QtChart import QChart, QChartView, QLineSeries, QVXYModelMapper, QAbstractSeries
from PyQt5.QtChart import QAbstractAxis, QValueAxis
from PyQt5.QtCore import QAbstractTableModel, QModelIndex, QRect, QEvent, QThread
from PyQt5.QtWidgets import QGridLayout, QHeaderView, QTableView, QWidget, QAction, QLineEdit, QLabel, QMessageBox, \
    QComboBox, QAbstractItemView
from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot, QPointF, Qt, QLineF, QPoint, pyqtBoundSignal, QStringListModel
from PyQt5.QtWidgets import QWidget, QApplication, QDialog, QMainWindow, QMenu, QVBoxLayout, QFileDialog, QPushButton
from PyQt5.QtGui import QPainter, QPolygonF, QActionEvent, QBrush, QPen, QColor, QPaintEvent, QMouseEvent
from PyQt5.QtGui import QFont, QKeyEvent, QFocusEvent, QIcon, QResizeEvent, QCloseEvent

from time import sleep
import numpy as np
import pandas as pd
from fluid import *
from core import *
from scipy import stats
import pickle as pkl
# from shutil import copyfile
from enum import Enum, auto
import sys

import produced_fluids


__version__ = '0.1.0'


def series_to_polyline(xdata, ydata):
    """Convert series data to QPolygon(F) polyline

    This code is derived from PythonQwt's function named 
    `qwt.plot_curve.series_to_polyline`"""

    if isinstance(xdata, np.ndarray):
        xdata = xdata.tolist()
    if isinstance(ydata, np.ndarray):
        ydata = ydata.tolist()

    size = len(xdata)
    polyline = QPolygonF(size)
    pointer = polyline.data()
    dtype, tinfo = np.float, np.finfo  # integers: = np.int, np.iinfo
    pointer.setsize(2*polyline.size()*tinfo(dtype).dtype.itemsize)
    memory = np.frombuffer(pointer, dtype)
    memory[:(size-1)*2+1:2] = xdata
    memory[1:(size-1)*2+2:2] = ydata
    return polyline


class CoreFloodSaveStructure:

    def __init__(self, cf):
        self.name = cf.name
        self.temperature = cf.temperature
        self.period = cf.period
        self.flow_rate = cf.flow_rate
        self.time_data = cf.time_data
        self.flow_data = cf.flow_data
        self.p_data = cf.p_data
        self.plateau_x = cf.plateau_x
        self.plateau_whole = cf.plateau_whole
        self.plateau_sec1 = cf.plateau_sec1
        self.plateau_sec2 = cf.plateau_sec2
        self.plateau_sec3 = cf.plateau_sec3
        self.plateau_sec4 = cf.plateau_sec4
        self.core = cf.core
        self.fluid = cf.fluid


class CoreFlood:

    # pd_change = pyqtSignal(str, bool, int)

    def __init__(self, name: str, my_fluid, my_core, my_experiment):
        super(CoreFlood, self).__init__()
        self.name = name
        self.notes = ''
        self.observations = ''
        self.source_file = ''
        self.fluid = my_fluid
        self.core = my_core
        self.experiment = my_experiment
        self.cut_num = 0
        self.flood_view = None
        self.temperature = 23.
        self.back_pressure = 0.
        self.planned_volume = 500.
        self.effluent = produced_fluids.ProducedFluids(self, True, True)
        # self.curves = \
        #     [QLineSeries(), QLineSeries(), QLineSeries(), QLineSeries(), QLineSeries()]
        self.period = 1.0
        self.flow_rate = {0.: 1.}  # Python "dictionary" of flow rate breaks and flow rates {break [min]: rate [mL/min]}
        self.time_data = np.ndarray([])  # raw time data [min]
        self.flow_data = np.ndarray([])  # volume data [mL]
        self.p_data = [np.ndarray([]), np.ndarray([]), np.ndarray([]), np.ndarray([]), np.ndarray([])]
        self.rf_data = [np.ndarray([]), np.ndarray([]), np.ndarray([]), np.ndarray([]), np.ndarray([])]
        self.plateau_x = [[], []]  # where plateau ends are stored [min]
        self.plateau_whole = []  # actual plateau averages [psi]
        self.plateau_sec1 = []
        self.plateau_sec2 = []
        self.plateau_sec3 = []
        self.plateau_sec4 = []
        self.plateau_whole_rf = []  # actual plateau rf averages
        self.plateau_sec1_rf = []
        self.plateau_sec2_rf = []
        self.plateau_sec3_rf = []
        self.plateau_sec4_rf = []

        self.pore_volume = np.nan
        self.tracer_object = None
        self.use_as_ref = False
        self.permeability = np.zeros(5)
        self.permeability[:] = np.nan
        self.permeability_viscosity = np.nan
        self.krw = np.nan
        self.kro = np.nan
        self.offsets = np.zeros(5)
        self.manual_offsets = np.zeros(5)

        self.used_offsets = self.manual_offsets
        self.offset_source_flood = self
        self.offset_referencer_floods = []
        self.permeability_source_flood = self
        self.permeability_referencer_floods = []
        self.rf_source_flood = self
        self.rf_referencer_floods = []
        self.reverse_flow = False
        self.n_phases = 1

        self.excluded = [[], [], [], [], []]
        self.y_lims = [np.nan, np.nan]
        self.y_lims_rf = [np.nan, np.nan]

    def get_experiment_total_pore_volume(self):

        return self.core.my_bulk_volume() * self.experiment.petro_parameters['phi'][1]

    def get_fluids_list(self) -> list:

        return [self.fluid]

    def get_last_fluid(self) -> object:

        return self.fluid

    def get_fluid_types(self) -> list:

        return [self.fluid.get_fluid_type()]

    def get_plateau_fluids(self, excluded) -> list:

        plateau_vs = self.get_plateau_volumes(excluded)
        if not plateau_vs:
            return []

        plateau_fluids = [self.fluid for i in range(len(plateau_vs))]
        # print('plateau volumes, fluids:', plateau_vs, plateau_fluids)

        return plateau_fluids

    def get_plateau_floods(self, excluded) -> list:

        plateau_vs = self.get_plateau_volumes(excluded)
        if not plateau_vs:
            return []

        plateau_floods = [self for i in range(len(plateau_vs))]

        return plateau_floods

    def get_rates_and_plateaus(self, i: int):
        if i == 1:
            plateaus = self.plateau_sec1
        elif i == 2:
            plateaus = self.plateau_sec2
        elif i == 3:
            plateaus = self.plateau_sec3
        elif i == 4:
            plateaus = self.plateau_sec4
        else:
            plateaus = self.plateau_whole

        if len(plateaus) < 1:
            return

        break_points = list(self.flow_rate.keys())
        break_points.sort()

        p_flow_rates = []

        for plateau, t in zip(plateaus, self.plateau_x[0]):
            # append a starter/placeholder flow rate for each pressure plateau
            p_flow_rates.append(self.flow_rate[break_points[0]])
            for break_point in break_points:
                # scan through break points and find the flow rate at each break point
                flow_rate = self.flow_rate[break_point]
                if break_point < t:
                    # if the break point occurs before the start of the plateau (t), update as the new flow rate
                    # the last flow rate to meet this criterion will actually be correct and supersede the others
                    p_flow_rates[-1] = flow_rate

        return p_flow_rates, plateaus

    def get_plateau_volumes(self, excluded: list) -> list:

        plateau_ts = self.plateau_x[0]
        plateau_ts = list(np.delete(plateau_ts, excluded))

        if not plateau_ts:
            return []

        plateau_vs = []
        my_pv = self.get_experiment_total_pore_volume()

        for t in plateau_ts:
            plateau_vs.append(self.map_t_to_pv(t) * my_pv)

        return plateau_vs

    def get_plateau_shear_rates(self, experiment, excluded: list):

        p_flow_rates, _ = self.get_rates_and_plateaus(0)
        p_flow_rates = np.delete(p_flow_rates, excluded)
        p_flow_rates = list(p_flow_rates)

        # p_floods = self.get_plateau_floods(excluded)

        us = np.divide(np.divide(p_flow_rates, 60.), self.core.my_area())
        # sw = []
        # for p_flood in p_floods:
        #     sw.append(experiment.get_flood_saturation(p_flood))
        # sw = np.array(sw)
        sw = experiment.get_flood_saturation(self)
        c = experiment.petro_parameters['C'][1]
        phi = experiment.petro_parameters['phi'][1]
        k = experiment.petro_parameters['perm'][1] * 0.000000009869 / 1000.
        swr = experiment.petro_parameters['Swr'][1]
        sor = experiment.petro_parameters['Sor'][1]
        s = (sw - swr) / (1. - swr - sor)

        if isinstance(self.fluid.specific_fluid, BrineInjectionFluid):
            if not np.isnan(self.krw):
                krw = self.krw
            else:
                krw0 = experiment.petro_parameters['krw0'][1]
                nw = experiment.petro_parameters['nw'][1]
                krw = krw0 * np.power(s, nw)

            shear = c * np.divide(us, np.sqrt(0.5 * np.multiply(k * krw, phi * sw)))

            if self.fluid.specific_fluid.polymer_solution is not None:
                n = self.fluid.specific_fluid.polymer_solution.model['n']
                if not np.isnan(n):
                    shear = ((3. * n + 1.) / (4. * n)) ** (n / (n - 1.)) * shear

            print('aqueous plateau k, krw, phi, s, sw: {:.1f}, {:.3f}, {:.4f}, {:.3f}, {:.3f}'.format(k, krw, phi, s, sw))
            print('aqueous plateau shear rates:', shear)
        else:
            if not np.isnan(self.kro):
                kro = self.kro
            else:
                kro0 = experiment.petro_parameters['kro0'][1]
                no = experiment.petro_parameters['no'][1]
                kro = kro0 * np.power(1. - s, no)

            shear = c * np.divide(us, np.sqrt(0.5 * np.multiply(k * kro, phi * (1. - sw))))

        return shear

    def add_data(self, x_data, y_data, n, color=None, mode='t', new_data=True):
        self.p_data[n] = y_data.copy()
        self.set_curve_data(x_data, y_data, n, color, mode, new_data)

    def set_curve_data(self, x_data, y_data, n, color=None, mode='t', new_data=True):

        self.time_data = x_data
        self.flood_view.x_min = x_data[0]
        self.flood_view.x_max = x_data[-1]
        self.set_flow_data()
        pen = self.flood_view.curves[n].pen()
        if color is not None:
            pen.setColor(color)
        pen.setWidthF(.5)
        self.flood_view.curves[n].setPen(pen)
        if mode == 'PV':
            pv = self.get_experiment_total_pore_volume()
            self.flood_view.curves[n].replace(series_to_polyline(self.flow_data / pv, y_data - self.used_offsets[n]))
        elif mode == 'V':
            self.flood_view.curves[n].replace(series_to_polyline(self.flow_data, y_data - self.used_offsets[n]))
        else:
            self.flood_view.curves[n].replace(series_to_polyline(x_data, y_data - self.used_offsets[n]))

        if self.flood_view is not None:
            print('setting axes, mode: ' + mode, new_data)
            self.flood_view.set_axes(mode, new_data, n)

    def set_curve_data_rf(self, x_data, y_data, n, color=None, mode='t', new_data=True):

        self.time_data = x_data
        self.flood_view.x_min = x_data[0]
        self.flood_view.x_max = x_data[-1]
        self.set_flow_data()
        pen = self.flood_view.curves[n].pen()
        if color is not None:
            pen.setColor(color)
        pen.setWidthF(.5)
        self.flood_view.curves[n].setPen(pen)

        if mode == 'PV':
            pv = self.get_experiment_total_pore_volume()
            self.flood_view.curves[n].replace(series_to_polyline(self.flow_data / pv, y_data))
        elif mode == 'V':
            self.flood_view.curves[n].replace(series_to_polyline(self.flow_data, y_data))
        else:
            self.flood_view.curves[n].replace(series_to_polyline(x_data, y_data))

        if self.flood_view is not None:
            print('setting axes, mode: ' + mode, new_data)
            self.flood_view.set_axes(mode, new_data, n)

    def set_flow_data(self):
        self.flow_data = self.time_data.copy()
        keys = list(self.flow_rate.keys())
        keys.sort()

        for key in keys:
            value = self.flow_rate[key]
            b = self.time_data < key
            prev_data = self.flow_data[b]
            prev_t = self.time_data[b]
            if len(prev_data) > 0:
                cum = prev_data[-1]
                t_offset = prev_t[-1]
            else:
                cum = 0.
                t_offset = 0.

            b = self.time_data >= key
            self.flow_data[b] = (self.time_data[b] - t_offset) * value + cum

    def get_flow_rate_vs_time(self) -> np.array:

        ts = list(self.flow_rate.keys())
        fvt = self.time_data.copy()

        for t in ts:
            mask = self.time_data >= t
            fvt[mask] = self.flow_rate[t]

        return fvt

    def update_rf_data(self):
        """ This routine updates the R.F. data stored in the flood object. """

        sf = self.rf_source_flood

        if sf == self:
            return

        if sf.permeability[0] == np.nan:
            return

        if sf.permeability_viscosity == np.nan:
            return

        # The calculation below is done in CGS units.

        # Get permeability in cm^2 and viscosity in P.
        k = sf.permeability * 9.86923E-12
        mu = sf.permeability_viscosity * 0.01

        # Get flow rate vs. time in cm^3/s, core length in cm, and core area in cm^2
        fvt = self.get_flow_rate_vs_time() / 60.
        l = self.experiment.core.my_length()
        a = self.experiment.core.my_area()

        for i in range(5):
            if i == 0:
                li = l
            else:
                if 7.62 * i < l:
                    li = 7.62
                else:
                    li = l - 7.62 * (i - 1.)
            # Get pressure data in Ba.
            p_data_offset = 68947.6 * (self.p_data[i] - self.used_offsets[i])
            # Get running implied permeability to viscosity ratio.
            perm_over_visc = np.divide(li * fvt, a * p_data_offset)
            # Resistance factor is ratio of perm to viscosity of reference flood to running implied perm to viscosity.
            self.rf_data[i] = np.divide(k[i] / mu, perm_over_visc)

    def update_plateaus(self):

        self.plateau_whole = []
        self.plateau_sec1 = []
        self.plateau_sec2 = []
        self.plateau_sec3 = []
        self.plateau_sec4 = []
        self.plateau_whole_rf = []
        self.plateau_sec1_rf = []
        self.plateau_sec2_rf = []
        self.plateau_sec3_rf = []
        self.plateau_sec4_rf = []
        rf_handles = [self.plateau_whole_rf, self.plateau_sec1_rf, self.plateau_sec2_rf, self.plateau_sec3_rf,
                      self.plateau_sec4_rf]

        if self.plateau_x[0]:
            for x1, x2 in zip(self.plateau_x[0], self.plateau_x[1]):
                if x1 > x2:
                    x = x1
                    x1 = x2
                    x2 = x
                td = self.time_data
                b = (td >= x1) & (td <= x2)  # create a boolean mask to index pressure data within the marked plateau

                self.plateau_whole.append(np.average(self.p_data[0][b]))  # select p_data in mask and take average [psi]
                self.plateau_sec1.append(np.average(self.p_data[1][b]))
                self.plateau_sec2.append(np.average(self.p_data[2][b]))
                self.plateau_sec3.append(np.average(self.p_data[3][b]))
                self.plateau_sec4.append(np.average(self.p_data[4][b]))

                if self.rf_data[0].size == self.p_data[0].size:
                    for i, rf_handle in enumerate(rf_handles):
                        b = (td >= x1) & (td <= x2) & ~np.isinf(self.rf_data[i])
                        if b.size:
                            rf_handle.append(np.average(self.rf_data[i][b]))
                        else:
                            rf_handle.append(np.nan)

    def refresh_curves_from_stored(self, x_data: list, mode):

        if self.flood_view is None:
            return

        if self.flood_view.y_view == FloodYView.PRESSURE:
            y_data_n = self.p_data.copy()
            for i, offset in enumerate(self.used_offsets):
                y_data_n[i] = y_data_n[i] - offset
        else:
            y_data_n = self.rf_data.copy()

        i = -1

        for curve, y_data in zip(self.flood_view.curves, y_data_n):
            i += 1
            yd = y_data.copy()
            curve.replace(series_to_polyline(x_data, yd.tolist()))

            if self.flood_view is not None:
                if mode == FloodXView.VOLUME:
                    pv = self.get_experiment_total_pore_volume()
                    self.flood_view.x_min = self.map_pv_to_t(min(x_data) / pv)
                    self.flood_view.x_max = self.map_pv_to_t(max(x_data) / pv)
                    self.flood_view.set_axes('V', False, 0)
                elif mode == FloodXView.PV:
                    self.flood_view.x_min = self.map_pv_to_t(min(x_data))
                    self.flood_view.x_max = self.map_pv_to_t(max(x_data))
                    self.flood_view.set_axes('PV', False, 0)
                else:
                    self.flood_view.x_min = min(x_data)
                    self.flood_view.x_max = max(x_data)
                    self.flood_view.set_axes('t', False, 0)

    def set_x_data_view_time(self):
        t_data = list(self.time_data.tolist())
        self.refresh_curves_from_stored(t_data, FloodXView.TIME)

    def set_x_data_view_volume(self):
        v_data = list(self.flow_data.tolist())
        self.refresh_curves_from_stored(v_data, FloodXView.VOLUME)

    def set_x_data_view_PV(self):
        pv = self.get_experiment_total_pore_volume()
        v_data = self.flow_data / pv
        v_data = list(v_data.tolist())
        self.refresh_curves_from_stored(v_data, FloodXView.PV)

    def map_pv_to_t(self, pv):

        i = np.arange(0, self.time_data.size, 1)
        c_pv = self.get_experiment_total_pore_volume()
        v = pv*c_pv
        index = np.interp(v, self.flow_data, i)
        t = np.interp(index, i, self.time_data)

        return t

    def map_t_to_pv(self, t):

        i = np.arange(0, self.time_data.size, 1)
        index = np.interp(t, self.time_data, i)
        c_pv = self.get_experiment_total_pore_volume()
        pv = np.interp(index, i, self.flow_data / c_pv)

        return pv

    def construct_pickle(self):
        print(self)

    def load_from_pickle(self, pickle: CoreFloodSaveStructure):
        self.name = pickle.name
        self.temperature = pickle.temperature
        self.period = pickle.period
        self.flow_rate = pickle.flow_rate
        self.time_data = pickle.time_data
        self.flow_data = pickle.flow_data
        self.p_data = pickle.p_data
        self.plateau_x = pickle.plateau_x
        self.plateau_whole = pickle.plateau_whole
        self.plateau_sec1 = pickle.plateau_sec1
        self.plateau_sec2 = pickle.plateau_sec2
        self.plateau_sec3 = pickle.plateau_sec3
        self.plateau_sec4 = pickle.plateau_sec4
        self.core = pickle.core
        self.fluid = pickle.fluid


class MultiCoreFlood(CoreFlood):

    def __init__(self, name: str, ref_objs: list, *args):
        my_fluid = ref_objs[0].fluid
        my_core = ref_objs[0].core
        super(MultiCoreFlood, self).__init__(name, my_fluid, my_core, my_experiment=args[0])

        self.ref_floods = ref_objs
        self.set_base_properties_from_flood(self.ref_floods[0])
        self.reset_effluent_from_ref_floods()

        self.view_class = MultiCoreFloodView

    def set_base_properties_from_flood(self, f: CoreFlood):

        self.temperature = f.temperature
        self.back_pressure = f.back_pressure
        self.n_phases = f.n_phases

    def reset_effluent_from_ref_floods(self):

        self.effluent = produced_fluids.ProducedFluids(self, True, True)
        vols = deepcopy(self.effluent.volumes)
        revs = deepcopy(self.effluent.revisions)
        ms = deepcopy(self.effluent.measurements)

        for fld in self.ref_floods:
            sv = fld.effluent.volumes.shape
            if sv[1]:
                sm = fld.effluent.measurements.shape
                s = vols.shape
                vols = np.resize(vols, (s[0], s[1] + sv[1]))
                vols[:, s[1]:] = fld.effluent.volumes
                s = revs.shape
                revs = np.resize(revs, (s[0], s[1] + sv[1]))
                revs[:, s[1]:] = fld.effluent.revisions
                s = ms.shape
                ms = np.resize(ms, (s[0], s[1] + sm[1]))
                ms[:, s[1]:] = fld.effluent.measurements

        self.effluent.volumes = vols
        self.effluent.measurements = ms

    def get_fluids_list(self) -> list:

        fluids = []
        for fl in self.ref_floods:
            fluids.append(fl.fluid)

        return fluids

    def get_last_fluid(self) -> InjectionFluid:

        return self.ref_floods[-1].fluid

    def get_fluid_types(self) -> list:

        fluid_types = []
        for fl in self.ref_floods:
            fluid_types.append(fl.fluid.get_fluid_type())

        return fluid_types

    def get_plateau_fluids(self, excluded) -> list:

        plateau_vs = self.get_plateau_volumes(excluded)
        if not plateau_vs:
            return []

        plateau_fluids = []

        for plateau_v in plateau_vs:
            v = 0.
            plateau_fluids.append(self.ref_floods[0].fluid)
            for fld in self.ref_floods:

                if plateau_v > v:
                    plateau_fluids[-1] = fld.fluid
                else:
                    break

                v += fld.planned_volume

        return plateau_fluids

    def get_plateau_floods(self, excluded) -> list:

        plateau_vs = self.get_plateau_volumes(excluded)
        if not plateau_vs:
            return []

        plateau_floods = []

        for plateau_v in plateau_vs:
            v = 0.
            plateau_floods.append(self.ref_floods[0])
            for fld in self.ref_floods:

                if plateau_v > v:
                    plateau_floods[-1] = fld
                else:
                    break

                v += fld.planned_volume

        return plateau_floods


class MultiCoreFloodView(QDialog):

    def __init__(self, mcf: MultiCoreFlood=None, parent=None):

        if mcf is None:
            return

        super(MultiCoreFloodView, self).__init__(parent=parent)
        self.setWindowTitle(' ')
        self.mcf = mcf

        lyt = QVBoxLayout()
        self.setLayout(lyt)

        name_label = QLabel(parent=self, text='Name: ' + mcf.name)
        ref_flood_listbox = QListWidget()
        for fld in self.mcf.ref_floods:
            ref_flood_listbox.addItem(fld.name)

        lyt.addWidget(name_label)
        lyt.addWidget(ref_flood_listbox)

        self.show()


class MultiCoreFloodBuilder(QDialog):

    def __init__(self, parent=None):
        super(MultiCoreFloodBuilder, self).__init__(parent=parent)
        # print(args)
        self.setWindowTitle('New MultiCoreFlood...')
        self.experiment = parent.parent().parent().experiment

        lyt = QVBoxLayout()
        self.setLayout(lyt)

        name_widget = QWidget(parent=self)
        name_lyt = QHBoxLayout()
        name_widget.setLayout(name_lyt)

        name_label = QLabel(parent=name_widget, text='Name: ')
        self.name_edit = QLineEdit(parent=name_widget)
        name_lyt.addWidget(name_label)
        name_lyt.addWidget(self.name_edit)

        self.floods_listbox = QListWidget()
        for fl in self.experiment.floods:
            self.floods_listbox.addItems([fl.name])
        self.floods_listbox.setSelectionMode(QAbstractItemView.ExtendedSelection)

        self.submit_button = QPushButton(parent=self, text='Submit')
        self.submit_button.clicked.connect(self.submit)

        lyt.addWidget(name_widget)
        lyt.addWidget(self.floods_listbox)
        lyt.addWidget(self.submit_button)

        self.show()

    def submit(self):

        name = self.name_edit.text()
        if not name:
            QMessageBox(parent=self, text='No name given.').exec_()
            return

        model_ind = self.floods_listbox.selectedIndexes()
        if not model_ind:
            QMessageBox(parent=self, text='No floods selected.').exec_()
            return

        inds = [m_ind.row() for m_ind in model_ind]

        if len(inds) > 1:
            for i, ind in enumerate(inds[:-1]):
                if ind + 1 != inds[i + 1]:
                    QMessageBox(parent=self, text='Floods must be consecutive.').exec_()
                    return

        flds = self.experiment.floods[inds[0]:inds[-1] + 1]

        for fld in flds:
            mf = self.experiment.find_multiflood(fld)
            if mf is not None:
                QMessageBox(parent=self, text='Flood ' + fld.name + ' is in use by multiflood ' + mf.name + '.').exec_()
                return
            if fld.fluid.get_base_fluid_type() != flds[0].fluid.get_base_fluid_type():
                QMessageBox(parent=self, text='All floods must use the same injection fluid base type.').exec_()
                return

        args = [self.experiment]
        print(self.parent(), name, flds, *args)
        try:
            self.parent().submit_name_ref_objects(name, flds, *args)
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()
        # self.parent().submit_name_ref_objects(name, flds, *args)
        self.close()


class CoreFloodFormSaveStructure:

    def __init__(self, cff):
        self.file_path = cff.file_path
        self.x_max = cff.x_max
        self.x_min = cff.x_min
        self.y_max = cff.y_max
        self.y_min = cff.y_min
        self.v_min = cff.v_min
        self.v_max = cff.v_max
        self.mode = cff.mode
        self.add = cff.add
        self.remove = cff.remove
        self.select = cff.select
        self.flood = cff.flood  # CoreFloodSaveStructure(cff.flood)


class Worker(QObject):

    finished = pyqtSignal()

    def __init__(self, flood_form):
        QObject.__init__(self)
        self.flood_form = flood_form
        self.stop_flag = False
        self.limit = 10000

    def start_continuous_update(self):
        self.flood_form.continuous_update_button.setText('Stop')
        self.flood_form.continuous_update_button.clicked.disconnect()
        self.flood_form.continuous_update_button.clicked.connect(self.flood_form.stop_continuous_update)
        self.flood_form.update_button.setEnabled(False)
        self.continuous_update(0)

    def continuous_update(self, i: int):

        while self.flood_form.continuous_update_button.text() == 'Stop' and i < self.limit:
            i += 1
            self.flood_form.update_from_source()
            if i == self.limit:
                self.flood_form.continuous_update_button.setText('Continuous')
            else:
                j = 0
                while j < 300:
                    j += 1
                    sleep(0.1)
                    if self.flood_form.continuous_update_button.text() == 'Continuous':
                        j = 300

        self.finished.emit()


class CoreFloodChartView(QChartView):

    def __init__(self, chart: QChart, linked_widget):
        QChartView.__init__(self, chart)
        self.linked_widget = linked_widget

    def mouseReleaseEvent(self, *args, **kwargs):
        QChartView.mouseReleaseEvent(self, *args, **kwargs)
        self.linked_widget.mouseReleaseEvent(*args)


class FloodXView(Enum):

    TIME = auto()
    VOLUME = auto()
    PV = auto()


class FloodYView(Enum):

    PRESSURE = auto()
    RF = auto()


class CoreFloodForm(QMainWindow):

    view_changed = pyqtSignal()

    def __init__(self, my_fluid, my_core, parent=None, linked_widget=None, flood=None, icon=None):
        super(CoreFloodForm, self).__init__(parent=parent)
        self.linked_widget = linked_widget
        self.flood_icon = icon
        if flood is None:
            self.flood = CoreFlood('dummy_flood', my_fluid, my_core)
        else:
            self.flood = flood

        self.flood.flood_view = self

        self.curves = \
            [QLineSeries(), QLineSeries(), QLineSeries(), QLineSeries(), QLineSeries()]

        self.v_lines = []
        self.plateau_lines = [[], []]
        self.file_path = ""
        self.chart = QChart(flags=Qt.WindowFlags())

        self.chart.legend().hide()
        icon = QIcon('UEORS logo cropped.png')
        self.setWindowIcon(icon)
        sf = linked_widget.width() / 1000
        self.view = CoreFloodChartView(self.chart, self)
        self.x_axis = QValueAxis()
        self.y_axis = QValueAxis()
        self.x_max = self.x_axis.max()
        self.x_min = self.x_axis.min()
        self.y_max = self.y_axis.max()
        self.y_min = self.y_axis.min()
        self.v_min = self.y_min
        self.v_max = self.y_max
        self.chart.setAxisY(self.y_axis)
        self.chart.setAxisX(self.x_axis)

        i = 0
        for curve in self.curves:
            if i == 0:
                curve.setName('Whole')
            else:
                curve.setName('Sec ' + str(i))
            i += 1
            self.chart.addSeries(curve)
            curve.attachAxis(self.x_axis)
            curve.attachAxis(self.y_axis)
            curve.setUseOpenGL(True)

        font = QFont()
        font.setPointSize(11)
        title_font = QFont()
        title_font.setPointSize(12)
        title_font.setBold(True)
        self.x_axis.setLabelsFont(font)
        self.y_axis.setLabelsFont(font)
        self.x_axis.setTitleFont(title_font)
        self.y_axis.setTitleFont(title_font)
        self.x_axis.setTitleText('Time [min]')
        self.y_axis.setTitleText(chr(916) + u'P [psi]')
        self.view.setRenderHint(QPainter.Antialiasing)

        self.view.setRubberBand(QChartView.RectangleRubberBand)

        self.toolbar = self.addToolBar("Navigation")
        self.toolbar.addAction('&Zoom', self.zoom_mode)
        self.toolbar.addAction('&Flow', self.flow_mode)
        self.toolbar.addAction('&Plateau', self.plateau_mode)
        self.zoom_toolbar = self.addToolBar("subbar")
        self.zoom_toolbar.addAction('&set limits', self.set_axis_limits)
        self.flow_toolbar = self.addToolBar("subbar2")
        self.flow_toolbar.addAction('&add', self.add_mode)
        self.flow_toolbar.addAction('&select/delete', self.select_mode)
        self.add = False
        self.remove = False
        self.select = False
        self.set_sub_toolbar()

        c_widget = QWidget()
        lyt = QVBoxLayout()
        l_widget = QWidget()
        l_widget.setFixedSize(int(sf * 800), int(sf * 20))

        self.legend_labels = [LegendLabel(l_widget, self, 0), LegendLabel(l_widget, self, 1), LegendLabel(l_widget, self, 2),
                              LegendLabel(l_widget, self, 3), LegendLabel(l_widget, self, 4)]
        self.legend_labels[0].setText('W')
        self.legend_labels[0].move(int(sf * 350), 0)
        l_colors = ['red', 'blue', 'green', 'magenta']
        l_font = self.legend_labels[0].font()
        # l_font.setBold(True)
        l_font.setPointSize(11)
        self.legend_labels[0].setFont(l_font)
        for i in range(4):
            self.legend_labels[i+1].setText(str(i+1))
            self.legend_labels[i+1].move(int(sf * (355+20*(i+1))), 0)
            self.legend_labels[i+1].setStyleSheet('QLabel { color : ' + l_colors[i] + '; }')
            self.legend_labels[i+1].setFont(l_font)

        # if isinstance(self.flood, MultiCoreFlood):
        #     f = self.flood.ref_floods[-1]
        # else:
        #     f = self.flood
        sw = self.linked_widget.experiment.get_flood_saturation(self.flood)
        saturation_text = 'Sw = %.3f' % sw
        self.saturation_label = QLabel(parent=l_widget, text=saturation_text)
        self.saturation_label.move(int(sf * 700), 0)
        self.saturation_label.setFont(l_font)

        self.l_widget = l_widget

        self.continuous_update_thread = None
        self.continuous_update_worker = None

        b_widget = QWidget()
        b_widget.setFixedSize(int(sf * 800), int(sf * 50))
        self.update_button = QPushButton('Update')
        self.update_button.clicked.connect(self.update_from_source)
        self.continuous_update_button = QPushButton('Continuous')
        self.continuous_update_button.clicked.connect(self.continuous_update_from_source)
        b_lyt = QHBoxLayout()
        b_lyt.setContentsMargins(0, 0, 0, 0)
        b_widget.setLayout(b_lyt)
        b_lyt.addWidget(self.update_button)
        b_lyt.addWidget(self.continuous_update_button)
        self.b_widget = b_widget
        lyt.addWidget(l_widget)
        lyt.addWidget(self.view)
        lyt.addWidget(self.b_widget)
        lyt.setAlignment(l_widget, Qt.AlignCenter)
        lyt.setAlignment(b_widget, Qt.AlignCenter)
        c_widget.setLayout(lyt)
        self.setCentralWidget(c_widget)
        if self.flood.source_file:
            file_list = self.flood.source_file.split('/')
            self.setWindowTitle('Core Flood Form: ' + self.flood.name + ' (' + file_list[-1] + ')')
        else:
            self.setWindowTitle('Core Flood Form: ' + self.flood.name)
        self.file_menu = QMenu('&File', self)
        self.file_menu.addAction('&Open', self.fileOpen)
        self.file_menu.addAction('&Save', self.file_save)
        self.file_menu.addAction('&Load', self.file_load)
        self.menuBar().addMenu(self.file_menu)
        self.view_menu = QMenu('&View', self)
        self.x_view_menu = QMenu('&x-view...')
        self.view_menu.addMenu(self.x_view_menu)
        self.x_view_menu.addAction('&Time [min]', self.view_time)
        self.x_view_menu.addAction('&Volume [mL]', self.view_volume)
        self.x_view_menu.addAction('&PV', self.view_PV)
        self.y_view_menu = QMenu('&y-view...')
        self.y_view_menu.addAction('&Pressure', self.view_pressure)
        self.y_view_menu.addAction('&R.F.', self.view_rf)
        self.view_menu.addMenu(self.y_view_menu)
        self.menuBar().addMenu(self.view_menu)
        self.tools_menu = QMenu('&Tools', self)
        self.tools_menu.addAction('&Calc. Perm', self.calc_perm)
        self.tools_menu.addAction('&Update R.F.', self.flood.update_rf_data)
        self.tools_menu.addAction('&Infer Visc.', self.infer_viscosity)
        self.tools_menu.addAction('&Effluent...', self.effluent)
        self.tools_menu.addAction('&Man. Offsets...', self.offsets)
        self.expt_params_menu = QMenu('&Expt Params...', self)
        self.tools_menu.addMenu(self.expt_params_menu)
        self.expt_params_menu.addAction('Use for perm...', self.use_for_perm)
        self.expt_params_menu.addAction('Use for endpoint...', self.use_for_endpoint)
        self.flood_params_menu = QMenu('&Flood Params...', self)
        self.flood_params_menu.addAction('Set krw...', self.set_krw)
        self.flood_params_menu.addAction('Set kro...', self.set_kro)
        self.tools_menu.addMenu(self.flood_params_menu)
        self.menuBar().addMenu(self.tools_menu)
        self.x_view = 'Time'
        self.y_view = FloodYView.PRESSURE
        self.mode = 'Zoom'
        self.x_axis.setVisible(False)
        self.y_axis.setVisible(False)
        self.statusBar().showMessage('Zoom...')

        self.rate_edit = CFFRateEdit(self)
        self.rate_edit.move(int(sf * 200), int(sf * 70))
        self.rate_edit.resize(int(sf * 30), int(sf * 30))
        self.rate_edit_label = QLabel(self)
        self.rate_edit_label.setText("mL/min")
        self.rate_edit_label.move(int(sf * 235), int(sf * 70))
        self.rate_edit.setVisible(False)
        self.rate_edit_label.setVisible(False)

        self.resize(int(sf * 800), int(sf * 700))

        self.effluent_view = None
        self.offsets_view = None
        self.krw_view = None
        self.kro_view = None
        self.ax_limits_dlg = None
        self.perm_view = None

        if self.flood.p_data[0].size < 2:
            self.activate_ui(False)
        else:
            self.activate_ui(True)
            self.x_axis.setVisible(True)
            self.y_axis.setVisible(True)

        if self.flood.time_data.size > 1:
            print('offsets:', self.flood.used_offsets)
            print('first data:', self.flood.p_data[0][0], self.flood.p_data[1][0])
            self.flood.set_curve_data(self.flood.time_data, self.flood.p_data[0], 0, Qt.black)
            self.flood.set_curve_data(self.flood.time_data, self.flood.p_data[1], 1, Qt.red)
            self.flood.set_curve_data(self.flood.time_data, self.flood.p_data[2], 2, Qt.blue)
            self.flood.set_curve_data(self.flood.time_data, self.flood.p_data[3], 3, Qt.darkGreen)
            self.flood.set_curve_data(self.flood.time_data, self.flood.p_data[4], 4, Qt.magenta)

            self.zoom_mode()

            self.x_max = self.x_axis.max()
            self.x_min = self.x_axis.min()
            if not np.isnan(self.flood.y_lims[0]):
                self.y_axis.setMin(self.flood.y_lims[0])
                self.y_axis.setMax(self.flood.y_lims[1])
            self.y_max = self.y_axis.max()
            self.y_min = self.y_axis.min()
            self.flood.y_lims = [self.y_min, self.y_max]
            self.v_max = self.y_max
            self.v_min = self.y_min

            self.re_add_lines()
            self.plateau_mode()

    def init_properties(self):
        # self.flood.pd_change.connect(self.set_axes)
        self.flood.flood_view = self
        self.v_lines = []
        self.plateau_lines = [[], []]
        self.file_path = ''

        i = 0
        for curve in self.curves:
            if i == 0:
                curve.setName('Whole')
            else:
                curve.setName('Sec ' + str(i))
            i += 1
            self.chart.addSeries(curve)
            curve.attachAxis(self.x_axis)
            curve.attachAxis(self.y_axis)
            curve.setUseOpenGL(True)

    def activate_ui(self, visible=True):
        self.toolbar.setVisible(visible)
        self.zoom_toolbar.setVisible(visible)
        self.flow_toolbar.setVisible(visible)
        self.l_widget.setVisible(visible)

    def set_axis_limits(self):
        if self.ax_limits_dlg is None:
            self.ax_limits_dlg = AxisLimitsDlg(self)

    def set_krw(self):
        if self.krw_view is None:
            try:
                self.krw_view = RelPermDlg(self, RelPermMode.krw)
            except Exception as e:
                print('Error with RelPermDlg:', e)

    def set_kro(self):
        if self.kro_view is None:
            try:
                self.kro_view = RelPermDlg(self, RelPermMode.kro)
            except Exception as e:
                print('Error with RelPermDlg:', e)

    def calc_perm(self):
        """This routine executes when the Calc Perm action from the Tools menu is clicked. It performs some checks on
        the pressure plateaus and fluids to screen whether to proceed with launching a permeability tool and if so,
        which tool to launch."""

        # Return if a perm view is already open for this flood.
        if self.perm_view is not None:
            return

        # Check first to see if there are any pressure plateaus. If not, return.
        r, p = self.flood.get_rates_and_plateaus(0)
        if not p:
            msg = QMessageBox(parent=self, text='Pressure plateaus required for permeability.')
            msg.exec_()
            return

        # Query if there is a selected plateau and find unique flow rates. If there are multiple plateaus but only one
        # unique rate and no selected plateau, abort.
        sel = self.get_selected_plateau()
        unique_r = list(np.unique(r))
        if len(p) > 1 and sel == -1 and len(unique_r) == 1:
            if len(p) > 1:
                msg = QMessageBox(parent=self, text='Plateau must be selected.')
                msg.exec_()
                return
            else:
                sel = 0

        # Query the fluids that are used in each pressure plateau. If no plateau is selected and analysis is to be run
        # on all plateaus, require the injection fluid to be the same for all.
        plateau_fluids = self.flood.get_plateau_fluids([])
        if sel == -1 and not all([pf == plateau_fluids[0] for pf in plateau_fluids]):
            return

        # Get a handle to the polymer_solution object and query fluid type.
        ps = plateau_fluids[sel].specific_fluid.polymer_solution
        fluid_type = plateau_fluids[sel].get_fluid_type()
        err_msg = ''

        # Do some important triage based on fluid_type. Live oils and oils or polymer solutions without viscosity data
        # should be excluded from the permeability analysis. Polymer solutions must also have a non-Newtonian model.
        if fluid_type == FluidType.LIVE_OIL:
            err_msg = 'Cannot use live oil for permeability.'

        elif fluid_type == FluidType.OIL:
            if not ps.stored_averages:
                err_msg = 'Oil has no shear viscosity data.'

        elif fluid_type == FluidType.POLYMER_SOLUTION:
            poly = ps.primary_additive
            b = ps.base_fluid
            c = ps.concentration
            if not ps.stored_averages and np.isnan(poly.estimate_viscosity(b, 30., 10., c)):
                err_msg = 'No viscosity data available for polymer solution.'
            elif ps.stored_averages and np.isnan(ps.model['n']):
                err_msg = 'Polymer solution is Newtonian.'
            elif not ps.stored_averages and len(unique_r) > 1 and sel == -1:
                err_msg = 'C-factor fitting requires shear curve for injection fluid.'

        elif fluid_type == FluidType.UNRECOGNIZED:
            err_msg = 'Fluid not recognized.'

        # If there's a non-empty error message, display it and abort.
        if err_msg:
            msg = QMessageBox(parent=self, text=err_msg)
            msg.exec_()
            return

        if fluid_type == FluidType.BRINE or fluid_type == FluidType.OIL or len(p) < 2 or sel != -1:
            # For brines or oils or polymer solutions when one plateau is available or selected, call the ordinary
            # PermView.
            try:
                self.perm_view = PermView(self, sel)
            except Exception as e:
                msg = QMessageBox(parent=self, text=str(e))
                msg.exec_()
                return

        elif fluid_type == FluidType.POLYMER_SOLUTION:
            # For polymer solutions when more than one plateau is available, call PolymerPermView.
            try:
                self.perm_view = PolymerPermView(self)
            except Exception as e:
                msg = QMessageBox(parent=self, text=str(e))
                msg.exec_()
                return

    def get_rates_and_plateaus(self, i: int):

        r, p = self.flood.get_rates_and_plateaus(i)
        return r, p

    def get_selected_plateau(self):

        sel = -1

        for i, p_line in enumerate(self.plateau_lines[0]):
            pen = p_line.pen()
            if pen.color() == Qt.yellow:
                sel = i
                break

        return sel

    def infer_viscosity(self):

        sf = self.flood.get_last_fluid().specific_fluid

        if isinstance(sf, BrineInjectionFluid) or sf.polymer_solution is not None:
            try:
                msg = QMessageBox(parent=self, text='Error: Viscosity can only be inferred for a live oil.')
                msg.exec_()
            except Exception as e:
                print(e)
            return
        elif isinstance(sf, OilInjectionFluid) and sf.polymer_solution is None:
            if self.flood.permeability_source_flood is self.flood:
                msg = QMessageBox(parent=self, text='Error: Permeability must reference another flood.')
                msg.exec_()
                return
            else:
                perm = self.flood.permeability[0]
                rates, plateaus = self.get_rates_and_plateaus(0)
                p_fluid = self.flood.get_plateau_fluids([])[-1]

                if p_fluid.specific_fluid != sf:
                    QMessageBox(
                        parent=self,
                        text='Error: Last plateau must be associated with last injection fluid (live oil).').exec_()
                    return
                rate = rates[-1]
                dP = plateaus[-1] - self.flood.used_offsets[0]
                A = self.flood.core.my_area()
                L = self.flood.core.my_length()
                mu = perm * A * dP / (245. * rate * L)
                sf.set_inferred_viscosity(mu)
                msg = QMessageBox(parent=self, text='Live Oil viscosity = ' + str(np.round(mu, 1)) + ' cP')
                msg.exec_()

    def effluent(self):

        if self.effluent_view is None:
            try:
                self.effluent_view = produced_fluids.ProducedFluidsView(self.flood.effluent, self)
            except Exception as e:
                QMessageBox(parent=self, text=str(e)).exec_()

    def offsets(self):

        if not hasattr(self, 'offsets_view'):
            self.__dict__['offsets_view'] = OffsetTool(self)
        elif self.offsets_view is None:
            self.offsets_view = OffsetTool(self)

    def use_for_perm(self):

        if np.isnan(self.flood.permeability[0]):
            msg = QMessageBox(parent=self, text='Permeability not defined.')
            msg.exec_()
            return

        pps = self.linked_widget.experiment.petro_parameters

        pps['perm'][1] = self.flood.permeability[0]
        pps['perm'][2] = self.flood

        if pps['krw0'][2] is not None:
            pps['krw0'][1] = pps['krw0'][2].permeability[0] / self.flood.permeability[0]

        if pps['kro0'][2] is not None:
            pps['kro0'][1] = pps['kro0'][2].permeability[0] / self.flood.permeability[0]

        if self.flood.n_phases != 1:
            msg = QMessageBox(parent=self,
                              text='Warning: generally, you should not use two-phase data for absolute permeability.')
            msg.exec_()

        if self.linked_widget.petro_view is not None:
            self.linked_widget.petro_view.estimates_tab.set_values_and_sources()

    def use_for_endpoint(self):

        if np.isnan(self.flood.permeability[0]):
            msg = QMessageBox(parent=self, text='Permeability not defined.')
            msg.exec_()
            return

        if self.flood.n_phases == 1:
            msg = QMessageBox(parent=self, text='Cannot set endpoint relative permeability using single-phase data.')
            msg.exec_()
            return

        sat = self.linked_widget.experiment.get_flood_saturation(self.flood)
        my_fluid = self.flood.get_last_fluid()

        if isinstance(my_fluid.specific_fluid, BrineInjectionFluid):
            self.linked_widget.experiment.petro_parameters['krw0'][1] = \
                self.flood.permeability[0] / self.linked_widget.experiment.petro_parameters['perm'][1]
            self.linked_widget.experiment.petro_parameters['krw0'][2] = self.flood
            self.linked_widget.experiment.petro_parameters['Sor'][1] = 1. - sat
            self.linked_widget.experiment.petro_parameters['Sor'][2] = self.flood

        else:
            self.linked_widget.experiment.petro_parameters['kro0'][1] = \
                self.flood.permeability[0] / self.linked_widget.experiment.petro_parameters['perm'][1]
            self.linked_widget.experiment.petro_parameters['kro0'][2] = self.flood
            # self.linked_widget.experiment.petro_parameters['krw0'][2] = self.flood
            self.linked_widget.experiment.petro_parameters['Swr'][1] = sat
            self.linked_widget.experiment.petro_parameters['Swr'][2] = self.flood
            self.linked_widget.experiment.petro_parameters['SoI'][1] = 1. - sat
            self.linked_widget.experiment.petro_parameters['SoI'][2] = self.flood

        if self.linked_widget.petro_view is not None:
            self.linked_widget.petro_view.estimates_tab.set_values_and_sources()

    def toggle_curve_visibility(self, index=0):
        if self.curves[index].isVisible():
            self.curves[index].setVisible(False)
        else:
            self.curves[index].setVisible(True)

    def focusInEvent(self, a0: QFocusEvent):
        print('focused!')

    def changeEvent(self, a0: QEvent):
        if a0.type() == QEvent.ActivationChange:
            print('activation changed!')
            if self.isActiveWindow():
                print('activated!')
            else:
                self.check_incomplete_plateau()
                print('deactivated!')

    def flow_mode(self):
        self.check_incomplete_plateau()
        self.view_time()
        self.add = False
        self.remove = False
        self.select = True
        self.v_lines_visible(True)
        self.plateau_lines_visible(False)
        self.view.setRubberBand(QChartView.NoRubberBand)
        self.set_sub_toolbar('flow')
        self.mode = 'Flow'
        self.force_chart_refresh()

    def zoom_mode(self):
        self.check_incomplete_plateau()
        self.add = False
        self.remove = False
        self.select = False
        self.v_lines_visible(False)
        self.plateau_lines_visible(False)
        self.view.setRubberBand(QChartView.RectangleRubberBand)
        self.set_sub_toolbar()
        self.mode = 'Zoom'
        self.force_chart_refresh()

    def plateau_mode(self):
        self.check_incomplete_plateau()
        self.view_PV()
        # self.view_menu.exec(QPoint(1, 1))
        self.add = False
        self.remove = False
        self.select = True
        self.v_lines_visible(False)
        self.plateau_lines_visible(True)
        self.view.setRubberBand(QChartView.NoRubberBand)
        self.set_sub_toolbar('flow')
        self.mode = 'Plateau'
        s = self.size()
        self.resize(s.width() + 1, s.height())
        self.resize(s.width(), s.height())
        self.force_chart_refresh()

    def add_mode(self):
        self.check_incomplete_plateau()
        self.add = True
        self.remove = False
        self.select = False
        self.v_lines_gray()
        if self.mode == 'Plateau':
            self.plateau_lines_visible(False)
            self.v_lines_visible(False)

    def remove_mode(self):
        self.check_incomplete_plateau()
        self.add = False
        self.remove = True
        self.select = False
        self.v_lines_gray()

    def select_mode(self):
        self.check_incomplete_plateau()
        self.add = False
        self.remove = False
        self.select = True
        if self.mode == 'Plateau':
            self.plateau_lines_visible(True)
            self.v_lines_visible(False)

    def v_lines_visible(self, v):
        if not v:
            self.rate_edit.setVisible(v)
            self.rate_edit_label.setVisible(v)
        for line in self.v_lines:
            line.setVisible(v)

    def v_lines_gray(self):
        self.rate_edit.setVisible(False)
        self.rate_edit_label.setVisible(False)
        for line in self.v_lines:
            pen = line.pen()
            pen.setColor(Qt.gray)
            line.setPen(pen)
            line.setVisible(False)
            line.setVisible(True)

    def plateau_lines_visible(self, v):
        for line1, line2 in zip(self.plateau_lines[0], self.plateau_lines[1]):
            line1.setVisible(v)
            line2.setVisible(v)

    def set_sub_toolbar(self, bar_name='zoom'):
        if bar_name == 'flow':
            self.flow_toolbar.setVisible(True)
            self.zoom_toolbar.setVisible(False)
        else:
            self.flow_toolbar.setVisible(False)
            self.zoom_toolbar.setVisible(True)

    def force_chart_refresh(self):
        maxi = self.x_axis.max()
        self.x_axis.setMax(maxi + 1)
        self.x_axis.setMax(maxi)

    def update_plateau_line_x(self):

        if not self.plateau_lines[0]:
            return

        for i in range(len(self.plateau_lines[0])):
            x1 = self.flood.plateau_x[0][i]
            x1 = self.flood.map_t_to_pv(x1)
            x2 = self.flood.plateau_x[1][i]
            x2 = self.flood.map_t_to_pv(x2)
            self.plateau_lines[0][i].clear()
            self.plateau_lines[1][i].clear()
            self.plateau_lines[0][i].append(series_to_polyline([x1, x1], [self.y_min, self.y_max]))
            self.plateau_lines[1][i].append(series_to_polyline([x2, x2], [self.y_min, self.y_max]))

    def add_v_line(self, x=0.):
        self.v_lines.append(QLineSeries())
        self.chart.addSeries(self.v_lines[-1])
        self.v_lines[-1].attachAxis(self.x_axis)
        self.v_lines[-1].attachAxis(self.y_axis)
        self.v_lines[-1].setUseOpenGL(True)

        self.v_lines[-1].append(series_to_polyline([x, x], [self.v_min, self.v_max]))
        pen = self.v_lines[-1].pen()
        pen.setColor(Qt.gray)
        pen.setWidthF(2)
        self.v_lines[-1].setPen(pen)
        self.update_plateau_line_x()

    def add_plateau_line(self, x=0.):
        index = 0
        if not self.plateau_lines[0]:
            self.plateau_lines[0].append(QLineSeries())
        elif len(self.plateau_lines[0]) == len(self.plateau_lines[1]):
            self.plateau_lines[0][-1].setVisible(False)
            self.plateau_lines[1][-1].setVisible(False)
            self.plateau_lines[0].append(QLineSeries())
        else:
            self.plateau_lines[1].append(QLineSeries())
            index = 1

        self.chart.addSeries(self.plateau_lines[index][-1])
        self.plateau_lines[index][-1].attachAxis(self.x_axis)
        self.plateau_lines[index][-1].attachAxis(self.y_axis)
        self.plateau_lines[index][-1].setUseOpenGL(True)

        self.plateau_lines[index][-1].append(series_to_polyline([x, x], [self.y_min, self.y_max]))
        pen = self.plateau_lines[index][-1].pen()
        pen.setColor(Qt.red)
        pen.setWidthF(2)
        self.plateau_lines[index][-1].setPen(pen)

        if index == 1:
            x1 = self.plateau_lines[0][-1].pointsVector()[-1].x()
            x2 = self.plateau_lines[1][-1].pointsVector()[-1].x()
            tx1 = self.flood.map_pv_to_t(x1)
            tx2 = self.flood.map_pv_to_t(x2)
            self.flood.plateau_x[0].append(tx1)
            self.flood.plateau_x[1].append(tx2)
            self.flood.update_plateaus()

    def check_incomplete_plateau(self):
        if len(self.plateau_lines[0]) != len(self.plateau_lines[1]):
            self.chart.removeSeries(self.plateau_lines[0][-1])
            self.plateau_lines[0].pop()

    def remove_v_line(self, v_line=None):
        if v_line is None:
            if not self.v_lines:
                pass
            else:
                series = self.chart.series()
                if self.v_lines[-1] in series:
                    print('found it')
                self.chart.removeSeries(self.v_lines[-1])
                self.v_lines[-1].detachAxis(self.x_axis)
                self.v_lines[-1].detachAxis(self.y_axis)
                self.v_lines.pop()
                self.force_chart_refresh()
        else:
            self.chart.removeSeries(v_line)
            self.v_lines.remove(v_line)
            self.force_chart_refresh()

    def remove_plateau_lines(self, p_line1, p_line2):
        self.chart.removeSeries(p_line1)
        self.chart.removeSeries(p_line2)
        self.plateau_lines[0].remove(p_line1)
        self.plateau_lines[1].remove(p_line2)
        self.force_chart_refresh()

    def set_title(self, title):
        self.chart.setTitle(title)

    def reset_axes(self, x_data, y_data):
        if min(x_data) < self.x_axis.min():
            self.x_axis.setMin(min(x_data))
        if max(x_data) > self.x_axis.max():
            self.x_axis.setMax(max(x_data))
        if min(y_data) < self.y_axis.min():
            self.y_axis.setMin(min(y_data))
            self.v_min = min(y_data)
        if max(y_data) > self.y_axis.max():
            self.y_axis.setMax(max(y_data))
            self.v_max = max(y_data)
        self.update_v_and_plateau_lines()

    def update_v_and_plateau_lines(self):
        for line in self.v_lines:
            pv = line.pointsVector()
            line.clear()
            line.append(series_to_polyline([pv[0].x(), pv[1].x()], [self.v_min, self.v_max]))
        for line in self.plateau_lines[0]:
            pv = line.pointsVector()
            line.clear()
            line.append(series_to_polyline([pv[0].x(), pv[1].x()], [self.v_min, self.v_max]))
        for line in self.plateau_lines[1]:
            pv = line.pointsVector()
            line.clear()
            line.append(series_to_polyline([pv[0].x(), pv[1].x()], [self.v_min, self.v_max]))

    def set_axes(self, mode, new_data=False, n=0):

        if not new_data:
            if mode == 't':
                print('here t:', self.x_min, self.x_max)
                self.x_axis.setMin(self.x_min)
                self.x_axis.setMax(self.x_max)
                # self.y_axis.setMin(self.y_min)
                # self.y_axis.setMax(self.y_max)
            elif mode == 'PV':
                f = self.flood
                print('here', self.x_max, f.map_t_to_pv(self.x_max))
                self.x_axis.setMin(f.map_t_to_pv(self.x_min))
                self.x_axis.setMax(f.map_t_to_pv(self.x_max))
                # self.y_axis.setMin(f.map_t_to_pv(self.y_min))
                # self.y_axis.setMax(f.map_t_to_pv(self.y_max))
            elif mode == 'V':
                f = self.flood
                pv = f.get_experiment_total_pore_volume()
                # pv = f.core.my_estimated_pore_volume()
                self.x_axis.setMin(pv*f.map_t_to_pv(self.x_min))
                self.x_axis.setMax(pv*f.map_t_to_pv(self.x_max))
                # self.y_axis.setMin(pv*f.map_t_to_pv(self.y_min))
                # self.y_axis.setMax(pv*f.map_t_to_pv(self.y_max))
        else:
            t_data = self.flood.time_data
            if self.y_view == FloodYView.PRESSURE:
                y_data = self.flood.p_data[n] - self.flood.used_offsets[n]
            else:
                y_data = self.flood.rf_data[n]
                non_inf = np.where(np.abs(y_data) != np.inf)
                y_data = y_data[non_inf]
                print('y_data limits for rf:', [min(y_data.tolist()), max(y_data.tolist())])
            self.reset_axes(t_data.tolist(), y_data.tolist())

    def update_from_source(self):
        try:
            self.x_axis.setVisible(True)
            self.y_axis.setVisible(True)
            txt = pd.read_csv(self.flood.source_file, skiprows=12, sep='\t')
            header = pd.read_csv(self.flood.source_file, nrows=0, skiprows=2)
        except Exception as e:
            msg = QMessageBox(parent=self, text=str(e))
            msg.exec_()
            return

        self.check_incomplete_plateau()
        period = float(str.split(header.columns[0])[-1]) / 60.0
        pressure_data = txt
        dims = pressure_data.shape
        size = dims[0]
        x_data = period * np.arange(0., size, 1.0)

        if self.x_view == 'PV':
            mode = 'PV'
        elif self.x_view == 'Volume':
            mode = 'V'
        else:
            mode = 't'
        print(mode)
        try:
            self.flood.add_data(x_data, pressure_data.values[:, 0], 0, color=Qt.black, mode=mode, new_data=False)
            self.flood.add_data(x_data, pressure_data.values[:, 1], 1, color=Qt.red, mode=mode, new_data=False)
            self.flood.add_data(x_data, pressure_data.values[:, 2], 2, color=Qt.blue, mode=mode, new_data=False)
            self.flood.add_data(x_data, pressure_data.values[:, 3], 3, color=Qt.darkGreen, mode=mode, new_data=False)
            self.flood.add_data(x_data, pressure_data.values[:, 4], 4, color=Qt.magenta, mode=mode, new_data=False)
        except Exception as e:
            QMessageBox(parent=self, text=str(e)).exec_()
            return

    def create_worker_and_thread(self):
        self.continuous_update_thread = QThread()
        self.continuous_update_worker = Worker(self)
        self.continuous_update_worker.moveToThread(self.continuous_update_thread)
        self.continuous_update_thread.started.connect(self.continuous_update_worker.start_continuous_update)
        self.continuous_update_worker.finished.connect(self.continuous_update_thread.quit)
        self.continuous_update_worker.finished.connect(self.continuous_update_worker.deleteLater)
        self.continuous_update_thread.finished.connect(self.continuous_update_thread.deleteLater)
        self.continuous_update_worker.finished.connect(lambda: self.continuous_update_button.setEnabled(True))
        self.continuous_update_worker.finished.connect(lambda: self.update_button.setEnabled(True))

    @pyqtSlot()
    def continuous_update_from_source(self):
        self.create_worker_and_thread()
        self.continuous_update_thread.start()

    @pyqtSlot()
    def stop_continuous_update(self):
        self.continuous_update_button.setText('Continuous')
        self.continuous_update_button.clicked.disconnect()
        self.continuous_update_button.clicked.connect(self.continuous_update_from_source)
        self.continuous_update_button.setEnabled(False)

    def load_data(self, pd_path, reset=True):
        self.x_axis.setVisible(True)
        self.y_axis.setVisible(True)
        txt = pd.read_csv(pd_path, skiprows=12, sep='\t')
        header = pd.read_csv(pd_path, nrows=0, skiprows=2)
        period = float(str.split(header.columns[0])[-1]) / 60.0
        pressure_data = txt
        dims = pressure_data.shape
        size = dims[0]
        x_data = period * np.arange(0., size, 1.0)
        self.flood.p_data = [np.ndarray([]), np.ndarray([]), np.ndarray([]), np.ndarray([]), np.ndarray([])]
        self.flood.plateau_x = [[], []]
        self.flood.plateau_whole = []
        self.flood.plateau_sec1 = []
        self.flood.plateau_sec2 = []
        self.flood.plateau_sec3 = []
        self.flood.plateau_sec4 = []
        self.flood.period = period
        self.flood.flow_rate = {0.: 1.}
        self.x_axis.setMin(0)
        self.x_axis.setMax(0)
        self.y_axis.setMin(0)
        self.y_axis.setMax(0)
        self.flood.add_data(x_data, pressure_data.values[:, 0], 0, color=Qt.black)
        self.flood.add_data(x_data, pressure_data.values[:, 1], 1, color=Qt.red)
        self.flood.add_data(x_data, pressure_data.values[:, 2], 2, color=Qt.blue)
        self.flood.add_data(x_data, pressure_data.values[:, 3], 3, color=Qt.darkGreen)
        self.flood.add_data(x_data, pressure_data.values[:, 4], 4, color=Qt.magenta)

        self.file_path = pd_path
        self.x_max = self.x_axis.max()
        self.x_min = self.x_axis.min()
        self.y_max = self.y_axis.max()
        self.y_min = self.y_axis.min()
        self.v_max = self.y_max
        self.v_min = self.y_min

        if self.v_lines:
            for line in self.v_lines:
                self.chart.removeSeries(line)

        if self.plateau_lines[0]:
            for line1, line2 in zip(self.plateau_lines[0], self.plateau_lines[1]):
                self.chart.removeSeries(line1)
                self.chart.removeSeries(line2)

        self.v_lines = []
        self.plateau_lines = [[], []]

        self.activate_ui(True)

        self.add_v_line()
        self.v_lines_visible(False)
        self.view_time()

    def fileQuit(self):
        self.close()

    def fileOpen(self):
        self.check_incomplete_plateau()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        # options |= QFileDialog.setDirectory('Q:')
        fileName, _ = QFileDialog.getOpenFileName(self, "Select pressure data",
                                                  self.file_path,
                                                  "Text Files (*.txt)", options=options)
        if fileName:
            try:
                self.load_data(fileName)
                self.flood.source_file = fileName
                file_list = self.flood.source_file.split('/')
                self.setWindowTitle('Core Flood Form: ' + self.flood.name + ' (' + file_list[-1] + ')')
            except Exception as e:
                msg = QMessageBox(parent=self, text=str(e))
                msg.exec_()

    def file_save(self):
        pickle = CoreFloodFormSaveStructure(self)
        i = self.file_path.rfind('/')
        j = self.file_path.rfind('.')
        file_name = self.file_path[i+1:j]
        pickle_out = open(file_name + '.flood', 'wb')
        pkl.dump(pickle, pickle_out)
        pickle_out.close()

    def file_load(self):

        dlg = QFileDialog()
        options = dlg.options()
        options |= dlg.DontUseNativeDialog
        dlg.setOptions(options)
        file_name, _ = dlg.getOpenFileName(self, "Select file", "",
                                           "Pickle Files (*.flood)")
        if not file_name:
            return

        pickle_in = open(file_name, 'rb')
        pickle = pkl.load(pickle_in)
        pickle_in.close()

        self.chart.removeAllSeries()
        core = Core('', [1.], [1.], [1.], 'sandstone')
        core.load_from_pickle(pickle.flood.core)
        self.flood = CoreFlood(pickle.flood.name, pickle.flood.fluid, core, pickle.flood.experiment)
        self.flood.flood_view = self
        self.init_properties()
        self.load_data(pickle.file_path)
        self.flood.load_from_pickle(pickle.flood)

        self.re_add_lines()

    def re_add_lines(self):

        self.flow_mode()
        for t in self.flood.flow_rate.keys():
            self.add_v_line(t)

        self.plateau_mode()
        plateau_x = self.flood.plateau_x
        self.flood.plateau_x = [[], []]
        for t1, t2 in zip(plateau_x[0], plateau_x[1]):
            v1 = self.flood.map_t_to_pv(t1)
            self.add_plateau_line(v1)
            v2 = self.flood.map_t_to_pv(t2)
            self.add_plateau_line(v2)

        self.zoom_mode()

    def construct_pickle(self):
        print(self)

    def closeEvent(self, ce: QCloseEvent):

        if self.krw_view is not None:
            self.krw_view.close()

        if self.ax_limits_dlg is not None:
            self.ax_limits_dlg.close()

        if self.flood is not None:
            self.flood.flood_view = None
        # self.flood_icon.cff = None

        try:
            print(self.flood_icon)
            exists = True
        except AttributeError:
            print('flood icon does not exist')
            exists = False

        try:
            if isinstance(self.flood, MultiCoreFlood):
                for i, fld in enumerate(self.flood.experiment.floods):
                    if fld in self.flood.ref_floods:
                        self.linked_widget.flood_forms[i] = None
                        self.linked_widget.flood_icons[i].cff = None

            else:
                if self.flood_icon is not None:
                    self.flood_icon.cff = None

                if self.linked_widget is not None:
                    ind = self.linked_widget.flood_icons.index(self.flood_icon)
                    self.linked_widget.flood_forms[ind] = None

        except Exception as e:
            print(e)

        # self.chart.removeAllSeries()
        if self.continuous_update_button.text() == 'Stop':
            self.continuous_update_button.setText('Continuous')
            sleep(0.2)

        self.destroy()

    def view_time(self):

        if self.flood.time_data.size == 1:
            return

        self.check_incomplete_plateau()
        self.x_axis.setMax(0.)
        # self.y_axis.setMin(self.y_min)
        # self.y_axis.setMax(self.y_max)
        self.flood.set_x_data_view_time()
        self.zoom_mode()
        self.x_view = 'Time'
        self.x_axis.setTitleText('Time [min]')
        self.view_changed.emit()
        self.v_lines_visible(False)

    def view_volume(self):

        if self.flood.time_data.size == 1:
            return

        self.check_incomplete_plateau()
        self.x_axis.setMax(0.)
        # self.y_axis.setMin(self.y_min)
        # self.y_axis.setMax(self.y_max)
        self.flood.set_x_data_view_volume()
        self.zoom_mode()
        self.v_lines_gray()
        self.x_view = 'Volume'
        self.x_axis.setTitleText('Volume [mL]')
        self.view_changed.emit()
        self.v_lines_visible(False)

    def view_PV(self):

        if self.flood.time_data.size == 1:
            return

        self.check_incomplete_plateau()
        self.x_axis.setMax(0.)
        # self.y_axis.setMin(self.y_min)
        # self.y_axis.setMax(self.y_max)
        self.flood.set_x_data_view_PV()
        self.zoom_mode()
        self.x_view = 'PV'
        self.v_lines_gray()
        self.x_axis.setTitleText('PV')
        self.view_changed.emit()
        self.v_lines_visible(False)

    def view_pressure(self):

        if self.y_view == FloodYView.PRESSURE:
            return

        if self.x_view == 'PV':
            mode = 'PV'
        elif self.x_view == 'Volume':
            mode = 'V'
        else:
            mode = 't'

        fl = self.flood

        for i in range(5):
            td = fl.time_data.copy()
            pd = fl.p_data[i].copy()
            fl.set_curve_data(td, pd, mode=mode, n=i, new_data=False)

        self.y_axis.setMin(self.flood.y_lims[0])
        self.y_axis.setMax(self.flood.y_lims[1])
        self.y_min = self.y_axis.min()
        self.y_max = self.y_axis.max()
        self.v_min = self.y_axis.min()
        self.v_max = self.y_axis.max()

        self.y_view = FloodYView.PRESSURE
        self.y_axis.setTitleText(chr(916) + u'P [psi]')

        self.update_v_and_plateau_lines()

    def view_rf(self):

        if self.y_view == FloodYView.RF:
            return

        fl = self.flood

        if not fl.p_data[0].size:
            QMessageBox(parent=self, text='No data available.').exec_()
            return

        fl_rf = fl.rf_source_flood

        if fl_rf == fl:
            QMessageBox(parent=self, text='Cannot calculate R.F. from self.').exec_()
            return

        if np.isnan(fl_rf.permeability[0]):
            QMessageBox(parent=self, text='Reference flood permeability not defined.').exec_()
            return

        if np.isnan(fl_rf.permeability_viscosity):
            QMessageBox(parent=self, text='Reference flood permeability viscosity not defined.').exec_()
            return

        rf_updated = False
        if fl.rf_data[0].size != fl.p_data[0].size:
            rf_updated = True
            self.flood.update_rf_data()

        if self.x_view == 'PV':
            mode = 'PV'
        elif self.x_view == 'Volume':
            mode = 'V'
        else:
            mode = 't'

        for i in range(5):
            td = fl.time_data.copy()
            rfd = fl.rf_data[i].copy()
            fl.set_curve_data_rf(td, rfd, i, mode=mode, new_data=False)

        self.y_view = FloodYView.RF
        self.y_axis.setTitleText('R.F.')

        if rf_updated:
            for i in range(5):
                self.set_axes(mode=mode, new_data=True, n=i)

            self.flood.y_lims_rf = [self.y_axis.min(), self.y_axis.max()]
            self.set_axes(mode=mode, new_data=False, n=0)
        else:
            self.y_axis.setMin(self.flood.y_lims_rf[0])
            self.y_axis.setMax(self.flood.y_lims_rf[1])

        self.y_min = self.y_axis.min()
        self.y_max = self.y_axis.max()
        self.v_min = self.y_axis.min()
        self.v_max = self.y_axis.max()

        self.update_v_and_plateau_lines()

        self.flood.update_plateaus()

    def mousePressEvent(self, a0: QMouseEvent):
        if a0.button() == 1:
            point = QPoint(a0.x()-10, a0.y()-90)
            scene_point = self.view.mapToScene(point)
            chart_point = self.chart.mapFromScene(scene_point)
            value_point = self.chart.mapToValue(chart_point)

            if value_point.x() > self.x_axis.max() or value_point.x() < self.x_axis.min():
                return

            if value_point.y() > self.y_axis.max() or value_point.y() < self.y_axis.min():
                return

            if self.add:
                if self.mode == 'Flow':
                    self.add_v_line(value_point.x())
                    self.flood.flow_rate[value_point.x()] = 1.
                    self.flood.set_flow_data()
                elif self.mode == 'Plateau':
                    self.add_plateau_line(value_point.x())

            elif self.remove:
                v_line = self.find_nearest_v_line(value_point.x())
                if v_line is not None:
                    x = v_line.pointsVector()[-1].x()
                    if x != 0:
                        if x in self.flood.flow_rate:
                            self.flood.flow_rate.pop(x)

                        self.remove_v_line(v_line)
            elif self.select:
                if self.mode == 'Flow':
                    if self.v_lines:
                        v_line = self.find_nearest_v_line(value_point.x())
                        for line in self.v_lines:
                            pen = line.pen()
                            if line == v_line:
                                p = line.pointsVector()[-1]
                                pp = self.chart.mapToPosition(p)
                                self.rate_edit.move(pp.x() + 5, 100)
                                self.rate_edit_label.move(pp.x() + 40, 100)
                                if pen.color() == Qt.yellow:
                                    pen.setColor(Qt.gray)
                                    self.rate_edit.setVisible(False)
                                    self.rate_edit_label.setVisible(False)
                                elif pen.color() == Qt.gray:
                                    pen.setColor(Qt.yellow)
                                    self.rate_edit.setVisible(True)
                                    self.rate_edit_label.setVisible(True)
                                    self.rate_edit.setText(str(self.flood.flow_rate[p.x()]))
                            else:
                                pen.setColor(Qt.gray)

                            line.setPen(pen)
                            line.setVisible(False)
                            line.setVisible(True)

                elif self.mode == 'Plateau':
                    if self.plateau_lines[0]:
                        p_line1, p_line2 = self.find_nearest_plateau_line(value_point.x())
                        self.statusBar().showMessage('')
                        i = -1
                        p_line_is = range(len(self.plateau_lines[0]))
                        for line1, line2 in zip(self.plateau_lines[0], self.plateau_lines[1]):
                            i += 1
                            pen = p_line1.pen()
                            if line1 == p_line1:
                                if pen.color() == Qt.yellow:
                                    pen.setColor(Qt.red)
                                    # plateau_flood = self.flood.get_plateau_floods([])[-1]
                                    # sw = self.linked_widget.experiment.get_flood_saturation(self.flood)

                                    # saturation_text = 'Sw = %.3f' % sw
                                    # self.saturation_label.setText(saturation_text)

                                elif pen.color() == Qt.red:
                                    pen.setColor(Qt.yellow)
                                    index = self.plateau_lines[0].index(line1)
                                    fl = self.flood

                                    excluded = np.setdiff1d(p_line_is, [i])
                                    shear = 0.
                                    try:
                                        shear = fl.get_plateau_shear_rates(self.linked_widget.experiment, excluded)[0]
                                    except Exception as e:
                                        msg = QMessageBox(parent=self, text=str(e))
                                        msg.exec_()

                                    mu = np.nan
                                    plateau_fluid = fl.get_plateau_fluids(excluded)[0]
                                    ps = plateau_fluid.specific_fluid.polymer_solution
                                    fluid_type = plateau_fluid.get_fluid_type()
                                    get_v = False

                                    if fluid_type == FluidType.BRINE or fluid_type == FluidType.LIVE_OIL:
                                        get_v = True
                                    elif fluid_type == FluidType.OIL and ps.stored_averages:
                                        get_v = True
                                    elif fluid_type == FluidType.POLYMER_SOLUTION and ps.stored_averages:
                                        get_v = True

                                    if get_v:
                                        try:
                                            print('temp, shear', fl.temperature, shear)
                                            mu = plateau_fluid.get_viscosity(fl.temperature, shear)
                                            print('mu:', mu)
                                        except Exception as e:
                                            msg = QMessageBox(parent=self, text=str(e))
                                            msg.exec_()

                                    elif fluid_type == FluidType.POLYMER_SOLUTION:
                                        poly = ps.primary_additive
                                        bf = ps.base_fluid
                                        conc = ps.concentration
                                        try:
                                            mu = poly.estimate_viscosity(bf, fl.temperature, shear, conc)
                                            if np.isnan(mu):
                                                mu = bf.get_viscosity(fl.temperature)
                                        except Exception as e:
                                            msg = QMessageBox(parent=self, text=str(e))
                                            msg.exec_()

                                    num_str = ''
                                    num_str += chr(947) + u' = %.1f s-1' % shear
                                    if not np.isnan(mu):
                                        num_str += u', ' + chr(956) + u' = %.1f cP' % mu
                                    else:
                                        num_str += u', ' + chr(956) + u' = NaN cP'

                                    if self.y_view == FloodYView.PRESSURE:
                                        num_str += u', whole = %.2f psi' % (fl.plateau_whole[index] - fl.used_offsets[0])
                                        num_str += u', sec1 = %.2f psi' % (fl.plateau_sec1[index] - fl.used_offsets[1])
                                        num_str += u', sec2 = %.2f psi' % (fl.plateau_sec2[index] - fl.used_offsets[2])
                                        num_str += u', sec3 = %.2f psi' % (fl.plateau_sec3[index] - fl.used_offsets[3])
                                        num_str += u', sec4 = %.2f psi' % (fl.plateau_sec4[index] - fl.used_offsets[4])
                                    else:
                                        num_str += u', whole = %.1f' % fl.plateau_whole_rf[index]
                                        num_str += u', sec1 = %.1f' % fl.plateau_sec1_rf[index]
                                        num_str += u', sec2 = %.1f' % fl.plateau_sec2_rf[index]
                                        num_str += u', sec3 = %.1f' % fl.plateau_sec3_rf[index]
                                        num_str += u', sec4 = %.1f' % fl.plateau_sec4_rf[index]

                                    self.statusBar().showMessage(num_str)

                                    # plateau_flood = fl.get_plateau_floods(excluded)[0]
                                    # sw = self.linked_widget.experiment.get_flood_saturation(plateau_flood)

                                    # saturation_text = 'Sw = %.3f' % sw
                                    # self.saturation_label.setText(saturation_text)
                            else:
                                pen.setColor(Qt.red)

                            line1.setPen(pen)
                            line1.setVisible(False)
                            line1.setVisible(True)
                            line2.setPen(pen)
                            line2.setVisible(False)
                            line2.setVisible(True)

    def mouseReleaseEvent(self, a0: QMouseEvent):
        if self.mode == 'Zoom' and a0.button() == 1:
            self.y_min = self.y_axis.min()
            self.y_max = self.y_axis.max()
            self.v_min = self.y_min
            self.v_max = self.y_max
            self.update_v_and_plateau_lines()

    def find_nearest_v_line(self, x):
        if not self.v_lines:
            return None
        else:
            nearest = self.v_lines[0]
            dx = nearest.pointsVector().pop()
            dx = abs(dx.x() - x)
            for v_line in self.v_lines:
                p = v_line.pointsVector()
                p = p.pop()
                if abs(x - p.x()) < dx:
                    nearest = v_line
                    dx = abs(x - p.x())
            return nearest

    def find_nearest_plateau_line(self, x):
        if not self.plateau_lines[0]:
            return None
        else:
            nearest1 = self.plateau_lines[0][0]
            dx1 = nearest1.pointsVector().pop()
            nearest2 = self.plateau_lines[1][0]
            dx2 = nearest2.pointsVector().pop()
            dx = abs(0.5*(dx1.x() + dx2.x()) - x)
            for plateau_line1, plateau_line2 in zip(self.plateau_lines[0], self.plateau_lines[1]):
                p1 = plateau_line1.pointsVector()
                p1 = p1.pop()
                p2 = plateau_line2.pointsVector()
                p2 = p2.pop()
                if abs(x - 0.5*(p1.x() + p2.x())) < dx:
                    nearest1 = plateau_line1
                    nearest2 = plateau_line2
                    dx = abs(x - 0.5*(p1.x() + p2.x()))
            return nearest1, nearest2

    def keyPressEvent(self, a0: QKeyEvent):

        if a0.key() == Qt.Key_Escape:

            self.x_min = self.flood.time_data[0]
            self.x_max = self.flood.time_data[-1]
            if self.x_view == 'PV':
                x_min = self.flood.map_t_to_pv(self.x_min)
                x_max = self.flood.map_t_to_pv(self.x_max)
            elif self.x_view == 'Volume':
                x_min = self.flood.map_t_to_pv(self.x_min) * self.flood.get_experiment_total_pore_volume()
                x_max = self.flood.map_t_to_pv(self.x_max) * self.flood.get_experiment_total_pore_volume()
            else:
                x_min = self.x_min
                x_max = self.x_max
            self.x_axis.setMax(x_max)
            self.x_axis.setMin(x_min)
            self.y_min = self.flood.y_lims[0]
            self.y_max = self.flood.y_lims[1]
            self.y_axis.setMax(self.y_max)
            self.y_axis.setMin(self.y_min)
            self.update_v_and_plateau_lines()

        elif a0.key() == Qt.Key_Delete:

            if self.mode == 'Plateau':
                if self.select:
                    if self.plateau_lines[0]:
                        for line1, line2 in zip(self.plateau_lines[0], self.plateau_lines[1]):
                            pen = line1.pen()
                            if pen.color() == Qt.yellow:
                                self.statusBar().showMessage('')
                                x1 = line1.pointsVector()[-1].x()
                                x2 = line2.pointsVector()[-1].x()
                                tx1 = self.flood.map_pv_to_t(x1)
                                tx2 = self.flood.map_pv_to_t(x2)
                                self.flood.plateau_x[0].remove(tx1)
                                self.flood.plateau_x[1].remove(tx2)
                                self.flood.update_plateaus()
                                self.remove_plateau_lines(line1, line2)

            elif self.mode == 'Flow':
                if self.select:
                    if self.v_lines:
                        for line in self.v_lines:
                            pen = line.pen()
                            if pen.color() == Qt.yellow:
                                x = line.pointsVector()[-1].x()

                                if x != 0:
                                    if x in self.flood.flow_rate:
                                        self.flood.flow_rate.pop(x)
                                    self.flood.set_flow_data()
                                    self.update_plateau_line_x()
                                    self.remove_v_line(line)
                                    self.rate_edit.setVisible(False)
                                    self.rate_edit_label.setVisible(False)

    def resizeEvent(self, a0: QResizeEvent):
        # super(CoreFloodForm, self).resizeEvent(a0=a0)
        self.correct_rate_edit_position()

    def correct_rate_edit_position(self):
        if self.v_lines:
            if self.v_lines[0].isVisible():
                for line in self.v_lines:
                    pen = line.pen()
                    if pen.color() == Qt.yellow:
                        pt = line.pointsVector()
                        x = pt.pop().x()
                        y = self.rate_edit.y()
                        new = self.chart.mapToPosition(QPoint(x, y))
                        self.rate_edit.move(new.x() + 5, 100)
                        self.rate_edit_label.move(new.x() + 40, 100)


class TempEdit(QLineEdit):

    def __init__(self, parent: QWidget, c_form: CoreFloodForm, x=730):
        super(TempEdit, self).__init__(parent)
        self.c_form = c_form
        self.setText(str(c_form.flood.temperature))
        self.setAlignment(Qt.AlignHCenter)
        self.resize(35, 20)
        self.move(x, 0)
        self.label = QLabel(parent)
        self.label.setText(chr(176) + u'C')
        self.label.resize(30, 20)
        self.label.move(x + 35, 0)

    def keyPressEvent(self, a0: QKeyEvent):
        super(TempEdit, self).keyPressEvent(a0)
        if a0.key() == Qt.Key_Enter:
            try:
                temp = float(self.text())
                self.c_form.flood.temperature = temp
            except ValueError:
                message = QMessageBox(self.c_form)
                message.setText('Invalid temperature.')
                message.show()


class CFFRateEdit(QLineEdit):

    def __init__(self, parent):
        super(CFFRateEdit, self).__init__(parent)
        self.editingFinished.connect(self.update_flow_rate)

    def update_flow_rate(self):

        p = self.parent()
        for line in p.v_lines:
            pen = line.pen()
            if pen.color() == Qt.yellow:
                x = line.pointsVector()[-1].x()
                p.flood.flow_rate[x] = float(self.text())
        p.flood.set_flow_data()
        p.update_plateau_line_x()


class LegendLabel(QLabel):

    label_clicked = pyqtSignal(int)

    def __init__(self, parent: QWidget, c_form: CoreFloodForm, index):
        super(LegendLabel, self).__init__(parent)
        self.label_clicked.connect(c_form.toggle_curve_visibility)

        self.index = index

    def mousePressEvent(self, ev: QMouseEvent):
        self.label_clicked.emit(self.index)


class AxisLimitsDlg(QDialog):

    def __init__(self, parent=None):
        super(AxisLimitsDlg, self).__init__(parent=parent)
        lyt = QVBoxLayout()
        self.edits = [QLineEdit(), QLineEdit(), QLineEdit(), QLineEdit()]
        lyt.addWidget(QLabel('x min'))
        lyt.addWidget(self.edits[0])
        lyt.addWidget(QLabel('x max'))
        lyt.addWidget(self.edits[1])
        lyt.addWidget(QLabel('y min'))
        lyt.addWidget(self.edits[2])
        lyt.addWidget(QLabel('y max'))
        lyt.addWidget(self.edits[3])
        self.btn = QPushButton('Enter')
        self.btn.clicked.connect(self.set_limits)
        # self.btn.connect(self.set_limits)
        lyt.addWidget(self.btn)
        self.setLayout(lyt)
        self.setFixedSize(150, 250)
        self.limits = [0, 1000, 0, 20]

        if parent is not None:
            self.correct_limits_for_view()
            parent.view_changed.connect(self.correct_limits_for_view)

        self.show()

    @pyqtSlot()
    def set_limits(self):

        limits = self.limits
        i = 0
        for edit in self.edits:
            try:
                limits[i] = float(edit.text())
                i += 1
            except ValueError:
                msg = QMessageBox(self)
                msg.setText('Invalid limit value.')
                msg.show()
                return

        p = self.parent()

        if p.x_view == 'Time':
            if limits[0] < p.flood.time_data[0]:
                limits[0] = p.flood.time_data[0]
                self.edits[0].setText(str(limits[0]))
            if limits[1] > p.flood.time_data[-1]:
                limits[1] = p.flood.time_data[-1]
                self.edits[1].setText(str(limits[1]))
            p.x_min = limits[0]
            p.x_max = limits[1]
        elif p.x_view == 'Volume':
            if limits[0] < p.flood.flow_data[0]:
                limits[0] = p.flood.flow_data[0]
                self.edits[0].setText(str(limits[0]))
            if limits[1] > p.flood.flow_data[-1]:
                limits[1] = p.flood.flow_data[-1]
                self.edits[1].setText(str(limits[1]))
            pv = p.flood.get_experiment_total_pore_volume()
            p.x_min = p.flood.map_pv_to_t(limits[0]/pv)
            p.x_max = p.flood.map_pv_to_t(limits[1]/pv)
        else:
            pv = p.flood.get_experiment_total_pore_volume()
            # pv = p.flood.core.my_estimated_pore_volume()
            if limits[0] < (p.flood.flow_data[0] / pv):
                limits[0] = p.flood.flow_data[0] / pv
                self.edits[0].setText(str(limits[0]))
            if limits[1] > (p.flood.flow_data[-1] / pv):
                limits[1] = p.flood.flow_data[-1] / pv
                self.edits[1].setText(str(limits[1]))
            p.x_min = p.flood.map_pv_to_t(limits[0])
            p.x_max = p.flood.map_pv_to_t(limits[1])

        p.x_axis.setMin(limits[0])
        p.x_axis.setMax(limits[1])
        p.y_axis.setMin(limits[2])
        p.y_axis.setMax(limits[3])

        p.y_min = limits[2]
        p.y_max = limits[3]
        p.v_min = limits[2]
        p.v_max = limits[3]
        if p.y_view == FloodYView.PRESSURE:
            p.flood.y_lims = [*limits[2:4]]
        else:
            p.flood.y_lims_rf = [*limits[2:4]]
        p.update_v_and_plateau_lines()

        p.correct_rate_edit_position()

        self.limits = limits

    def correct_limits_for_view(self):

        p = self.parent()
        self.limits = [p.x_min, p.x_max, p.y_min, p.y_max]

        if p.x_view == 'Time':
            pass
        elif p.x_view == 'Volume':
            pv = p.flood.get_experiment_total_pore_volume()
            # pv = p.flood.core.my_estimated_pore_volume()
            self.limits[0] = p.flood.map_t_to_pv(self.limits[0]) * pv
            self.limits[1] = p.flood.map_t_to_pv(self.limits[1]) * pv
        else:
            self.limits[0] = p.flood.map_t_to_pv(self.limits[0])
            self.limits[1] = p.flood.map_t_to_pv(self.limits[1])

        i = 0
        for edit in self.edits:
            edit.setText(str(self.limits[i]))
            i += 1

    def closeEvent(self, a0: QtGui.QCloseEvent):

        try:
            self.parent().ax_limits_dlg = None
        except Exception as e:
            print(e)


class RelPermMode(Enum):

    krw = auto()
    kro = auto()


class RelPermDlg(QDialog):

    def __init__(self, parent, mode: RelPermMode=RelPermMode.krw):
        super(RelPermDlg, self).__init__(parent=parent)
        self.mode = mode

        lyt = QVBoxLayout()

        fl = self.parent().flood
        if mode == RelPermMode.krw:
            perm_handle = fl.krw
            perm_text = 'krw = '
        else:
            perm_handle = fl.kro
            perm_text = 'kro = '

        self.perm_text = perm_text

        if not np.isnan(perm_handle):
            perm_text += str(perm_handle)

        self.perm_label = QLabel(perm_text)
        lyt.addWidget(self.perm_label)

        self.perm_edit = QLineEdit(parent=self)
        self.perm_edit.editingFinished.connect(self.check_entry)
        lyt.addWidget(self.perm_edit)

        self.set_button = QPushButton(parent=self, text='Set')
        self.set_button.clicked.connect(self.set_button_clicked)
        lyt.addWidget(self.set_button)
        self.clear_button = QPushButton(parent=self, text='Clear')
        self.clear_button.clicked.connect(self.clear_button_clicked)
        lyt.addWidget(self.clear_button)

        self.setLayout(lyt)
        print(self.size())
        self.show()

    def check_entry(self):

        txt = self.perm_edit.text()
        if not txt:
            return

        try:
            val = float(txt)
            if val < 0.:
                raise ValueError

        except ValueError:
            self.perm_edit.setText('')

    def set_button_clicked(self):

        self.check_entry()

        txt = self.perm_edit.text()
        if not txt:
            return

        val = float(txt)
        if self.mode == RelPermMode.krw:
            self.parent().flood.krw = val
        else:
            self.parent().flood.kro = val
        self.perm_label.setText(self.perm_text + str(val))

    def clear_button_clicked(self):

        if self.mode == RelPermMode.krw:
            self.parent().flood.krw = np.nan
        else:
            self.parent().flood.kro = np.nan

        self.perm_label.setText(self.perm_text)

    def closeEvent(self, a0: QtGui.QCloseEvent):
        print(a0.InputMethod)
        try:
            if self.mode == RelPermMode.krw:
                self.parent().krw_view = None
            else:
                self.parent().kro_view = None
        except Exception as e:
            print('RelPermDlg:', e)


class PermView(QDialog):

    def __init__(self, parent: CoreFloodForm, sel: int=-1):
        super(PermView, self).__init__(parent=parent)

        self.sel = sel

        if len(parent.flood.plateau_whole) == 0:
            del self
            return

        self.setWindowTitle('Permeability Tool')
        layout = QVBoxLayout()
        self.setLayout(layout)
        h_layout = QHBoxLayout()
        self.perm_text = QLabel(parent=self, text='Permeability: ')
        font = self.perm_text.font()
        font.setPointSize(12)
        self.perm_text.setFont(font)
        self.prev_button = QPushButton(parent=self, text='<')
        self.prev_button.setFixedWidth(20)
        self.prev_button.clicked.connect(self.prev_click)
        self.next_button = QPushButton(parent=self, text='>')
        self.next_button.setFixedWidth(20)
        self.next_button.clicked.connect(self.next_click)
        h_layout.addWidget(self.prev_button)
        h_layout.addWidget(self.perm_text)
        h_layout.addWidget(self.next_button)
        layout.addLayout(h_layout)
        self.i = 0
        self.color = Qt.black
        self.ex_color = Qt.gray

        if hasattr(parent.flood, 'excluded'):
            self.excluded = parent.flood.excluded
        else:
            self.excluded = [[], [], [], [], []]

        if len(parent.flood.plateau_whole) == 1 or sel != -1:

            self.setFixedSize(400, 50)
            self.chart = None

            self.calc_perm()

        else:

            self.setFixedSize(450, 450)

            self.chart = QChart(flags=Qt.WindowFlags())
            self.chart.legend().hide()
            self.view = QChartView(self.chart)

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
            self.x_axis.setTitleText('Q [mL/min]')
            self.y_axis.setTitleText(chr(916) + u'P [psi]')

            self.data_lines = []

            self.fit_line = QLineSeries()
            self.fit_line.setColor(self.color)
            self.fit_line.setUseOpenGL(False)
            self.chart.addSeries(self.fit_line)
            self.fit_line.attachAxis(self.x_axis)
            self.fit_line.attachAxis(self.y_axis)

            layout.addWidget(self.view)

            r, p = parent.flood.get_rates_and_plateaus(0)

            for rate, plateau in zip(r, p):
                data_line = QScatterSeries()
                data_line.setColor(self.color)
                data_line.setMarkerShape(QScatterSeries.MarkerShapeCircle)
                data_line.setMarkerSize(10.)
                data_line.setUseOpenGL(False)
                self.chart.addSeries(data_line)
                data_line.attachAxis(self.x_axis)
                data_line.attachAxis(self.y_axis)
                point = series_to_polyline([rate], [plateau])
                data_line.append(point)
                self.data_lines.append(data_line)

            self.x_axis.setMin(0.)
            self.x_axis.setMax(max(r) + min(r))
            self.fit_line_x = [0., max(r) + min(r)]
            self.y_axis.setMin(0.)
            self.y_axis.setMax(max(p) + 1.)

            self.calc_perm()

        for i in range(4):
            self.next_click()

        self.i = 1
        self.prev_click()

        self.show()

    def next_click(self):

        if self.i < 4:
            self.i += 1
            self.update_data()
            self.calc_perm()

    def prev_click(self):

        if self.i > 0:
            self.i -= 1
            self.update_data()
            self.calc_perm()

    def update_data(self):

        r, p = self.parent().flood.get_rates_and_plateaus(self.i)

        if self.sel != -1:
            r = [r[self.sel]]
            p = [p[self.sel]]

        if self.i == 1:
            self.color = Qt.red
            self.perm_text.setStyleSheet('color: red')
        elif self.i == 2:
            self.color = Qt.blue
            self.perm_text.setStyleSheet('color: blue')
        elif self.i == 3:
            self.color = Qt.darkGreen
            self.perm_text.setStyleSheet('color: green')
        elif self.i == 4:
            self.color = Qt.magenta
            self.perm_text.setStyleSheet('color: magenta')
        else:
            self.color = Qt.black
            self.perm_text.setStyleSheet('color: black')

        if len(p) < 2:
            return

        self.fit_line.setColor(self.color)

        j = -1
        for rate, plateau, line in zip(r, p, self.data_lines):
            j += 1
            line.clear()
            point = series_to_polyline([rate], [plateau])
            line.append(point)
            if j in self.excluded[self.i]:
                line.setColor(self.ex_color)
            else:
                line.setColor(self.color)

    def calc_perm(self):

        fv = self.parent()
        p_flow_rates, plateaus = fv.flood.get_rates_and_plateaus(self.i)

        if self.sel == -1:
            p_flow_rates = np.delete(np.array(p_flow_rates), self.excluded[self.i])
            plateaus = np.delete(np.array(plateaus), self.excluded[self.i])
            p_flow_rates = list(p_flow_rates.tolist())
            plateaus = list(plateaus.tolist())

        else:
            p_flow_rates = [p_flow_rates[self.sel]]
            plateaus = [plateaus[self.sel]]

        my_flood = self.parent().flood
        my_fluid = my_flood.get_plateau_fluids([])[self.sel]
        print(my_fluid.name)

        if len(p_flow_rates) > 1 and (isinstance(my_fluid.specific_fluid, OilInjectionFluid) or
                                              my_fluid.specific_fluid.polymer_solution is None):
            slope, intercept, r_value, p_value, std_err = stats.linregress(p_flow_rates, plateaus)
        elif len(p_flow_rates) == 1:
            intercept = self.parent().flood.used_offsets[self.i]
            slope = (plateaus[0] - intercept) / p_flow_rates[0]
        else:
            intercept = self.parent().flood.used_offsets[self.i]
            slope = np.divide(np.subtract(plateaus, intercept), p_flow_rates)

        T = fv.flood.temperature
        # T = 9. * T / 5. + 32.  # convert temperature to fahrenheit

        try:
            print(self.excluded)
            shear = fv.flood.get_plateau_shear_rates(fv.linked_widget.experiment, self.excluded[self.i])
            print(shear)
            # sw = fv.linked_widget.experiment.get_flood_saturation(my_flood)
            mu = my_fluid.get_viscosity(T, shear[0])
            # shear_msg = 'sw = {:.3f}, shear = {:.1f} s-1, mu = {:.1f} cP'.format(sw, shear[0], mu)
            # msg = QMessageBox(parent=self, text=shear_msg)
            # msg.exec_()
        except Exception as e:
            msg = QMessageBox(parent=self, text=str(e))
            msg.exec_()
            return

        c = fv.flood.core
        L = c.my_length()

        if self.i > 0:
            L -= (self.i - 1) * 7.62
            if L >= 7.62:
                L = 7.62

        if L < 0:
            self.prev_click()
            return

        A = c.my_area()

        if isinstance(slope, np.ndarray):
            perm = np.divide(245 * mu * L, A * slope)
            perm = float(perm[0])
        else:
            perm = 245. * mu * L / (A * slope)

        self.perm_text.setText('Permeability: ' + str(np.round(perm, 1)) + ' [md], Offset: ' +
                                str(np.round(intercept, 2)) + ' [psi]')

        fv.flood.permeability[self.i] = perm
        fv.flood.permeability_viscosity = mu
        # print('updating offsets in PermView.calc_perm().')
        fv.flood.offsets[self.i] = intercept

        if self.chart is not None:
            # print('updating chart.')
            self.y_axis.setMin(intercept)
            fit_line_y = [intercept, intercept + slope * self.fit_line_x[1]]
            self.y_axis.setMax(fit_line_y[1])
            fit_line_points = series_to_polyline(self.fit_line_x, fit_line_y)
            self.fit_line.clear()
            self.fit_line.append(fit_line_points)

    def mousePressEvent(self, a0: QtGui.QMouseEvent):

        if a0.button() != 1:
            return

        point = QPoint(a0.x(), a0.y())
        point = self.view.mapFromParent(point)
        point = self.chart.mapToValue(point)
        x = point.x()
        y = point.y()

        if x < self.x_axis.min() or x > self.x_axis.max() or y < self.y_axis.min() or y > self.y_axis.max():
            return

        try:
            r, p = self.parent().flood.get_rates_and_plateaus(self.i)
        except Exception as e:
            msg = QMessageBox(parent=self, text=str(e))
            msg.exec_()
            return

        index = utils.find_nearest(np.array(r), np.array(p), x, y)

        try:
            if self.data_lines[index].color() == self.color:
                if len(self.excluded[self.i]) < len(r) - 1:
                    self.data_lines[index].setColor(self.ex_color)
                    self.excluded[self.i].append(index)
            else:
                self.data_lines[index].setColor(self.color)
                self.excluded[self.i].remove(index)

        except Exception as e:
            msg = QMessageBox(parent=self, text=str(e))
            msg.exec_()
            return

        try:
            self.calc_perm()
            self.force_refresh()
        except Exception as e:
            msg = QMessageBox(parent=self, text=str(e))
            msg.exec_()
            return

    def force_refresh(self):

        w = self.width()
        self.setFixedWidth(w + 1)
        self.setFixedWidth(w)

        xm = self.x_axis.max()
        self.x_axis.setMax(xm + 1.)
        self.x_axis.setMax(xm)

    def closeEvent(self, a0: QtGui.QCloseEvent):

        try:
            fv = self.parent()
            fl = fv.flood
            for i in range(5):
                td = fl.time_data.copy()
                pd = fl.p_data[i].copy()
                fl.set_curve_data(td, pd, i)
            if hasattr(fl, 'excluded'):
                fl.excluded = self.excluded
            else:
                fl.__dict__['excluded'] = self.excluded

            p_flow_rates, _ = fv.flood.get_rates_and_plateaus(0)
            if len(p_flow_rates) > 1 and fl.offset_source_flood == fl:
                fl.used_offsets = fl.offsets

            if fv.mode == 'Zoom':
                fl.set_x_data_view_time()
                fv.zoom_mode()
            elif fv.mode == 'Plateau':
                fl.set_x_data_view_PV()
                fv.plateau_mode()
            else:
                fl.set_x_data_view_time()
                fv.flow_mode()

            fv.perm_view = None

        except Exception as e:
            print(e)


class PolymerPermView(PermView):

    def __init__(self, parent: CoreFloodForm):
        super(PolymerPermView, self).__init__(parent)

        my_flood = parent.flood
        self.mod_krw = False
        if not np.isnan(my_flood.krw):
            self.mod_krw = True

        self.params_widget = QWidget(parent=self)
        lyt = QHBoxLayout()
        self.params_widget.setLayout(lyt)

        self.c_text = QLabel(parent=self, text='C=')
        self.c_text.setFixedHeight(20)
        self.c_edit = QLineEdit(parent=self)
        self.c_edit.setFixedSize(25, 20)
        self.c_edit.setText(str(np.round(parent.linked_widget.experiment.petro_parameters['C'][1], 1)))
        self.c_edit.editingFinished.connect(self.update_params)
        self.set_c_push = QPushButton(parent=self, text='Set C')
        self.set_c_push.setFixedHeight(20)
        self.set_c_push.clicked.connect(self.set_c)

        self.krw0_text = QLabel(parent=self)
        self.krw0_text.setFixedHeight(20)
        self.krw0_edit = QLineEdit(parent=self)
        self.krw0_edit.setFixedSize(25, 20)
        if not self.mod_krw:
            krw0 = parent.linked_widget.experiment.petro_parameters['krw0'][1]
            krw_txt = 'krw0'
        else:
            krw0 = my_flood.krw
            krw_txt = 'krw'

        self.krw0_text.setText(krw_txt + '=')
        log_level = np.floor(np.log10(krw0))
        self.krw0_edit.setText(str(np.round(krw0 * 10. ** -log_level, 1)))
        self.krw0_logtext = QLabel(parent=self, text='x 10^')
        self.krw0_logtext.setFixedHeight(20)
        self.krw0_logedit = QLineEdit(parent=self)
        self.krw0_logedit.setFixedSize(25, 20)
        self.krw0_logedit.setText(str(np.round(log_level, 0)))
        self.set_krw0_push = QPushButton(parent=self, text='Set ' + krw_txt)
        self.set_krw0_push.setFixedHeight(20)
        self.set_krw0_push.clicked.connect(self.set_krw0)

        lyt.addWidget(self.c_text)
        lyt.addWidget(self.c_edit)
        lyt.addWidget(self.set_c_push)
        lyt.addWidget(self.krw0_text)
        lyt.addWidget(self.krw0_edit)
        lyt.addWidget(self.krw0_logtext)
        lyt.addWidget(self.krw0_logedit)
        lyt.addWidget(self.set_krw0_push)
        self.layout().addWidget(self.params_widget)

        self.fit_line.clear()
        self.chart.removeAllSeries()

        self.chart.removeAxis(self.x_axis)
        self.chart.removeAxis(self.y_axis)
        del self.x_axis
        del self.y_axis
        self.x_axis = QLogValueAxis()  # axes
        self.y_axis = QLogValueAxis()
        self.chart.setAxisX(self.x_axis)
        self.chart.setAxisY(self.y_axis)

        my_fluid = my_flood.get_plateau_fluids([])[-1]
        shear = np.logspace(-1., 3.)
        model = my_fluid.specific_fluid.polymer_solution.model
        visc = my_fluid.specific_fluid.polymer_solution.rheology_model(shear, model["eta_o"], model["eta_inf"],
                                                                  model["n"], model["l"])

        model_line = QLineSeries()
        model_line.setUseOpenGL(False)
        pen = model_line.pen()
        pen.setColor(Qt.black)
        pen.setWidthF(1.)
        model_line.setPen(pen)
        self.chart.addSeries(model_line)
        model_line.attachAxis(self.x_axis)
        model_line.attachAxis(self.y_axis)
        line = series_to_polyline(list(shear.tolist()), list(visc.tolist()))
        model_line.append(line)
        del self.fit_line
        self.fit_line = model_line
        self.x_axis.setMin(0.1)
        self.x_axis.setMax(1000.)
        self.x_axis.setTitleText('shear rate [1/s]')
        self.y_axis.setMin(1.0)
        self.y_axis.setMax(10. ** (np.log10(np.max(visc)) + 0.1))
        self.y_axis.setTitleText(chr(951) + u' [cP]')

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

        self.prev_button.setEnabled(False)
        self.next_button.setEnabled(False)

        try:
            self.update_data_lines()
            self.set_perm_text()
        except Exception as e:
            print(e)

    def calc_perm(self):
        print('PolymerPermView calc perm.', self)

    def set_perm_text(self):
        pps = self.parent().linked_widget.experiment.petro_parameters
        k = pps['perm'][1]
        if self.mod_krw:
            krw0 = self.parent().flood.krw
        else:
            krw0 = pps['krw0'][1]
        self.perm_text.setText('Permeability = ' + str(np.round(krw0 * k, 1)) + ' [md]')

    def update_params(self):

        params = self.parent().linked_widget.experiment.petro_parameters

        try:
            c = float(self.c_edit.text())
            assert c >= 0.
            params['C'][1] = np.round(c, 1)
            # self.update_data_lines()

        except Exception as e:
            msg = QMessageBox(parent=self, text=str(e))
            msg.exec_()
            return

        try:
            krw0 = float(self.krw0_edit.text())
            krw0_log = np.round(float(self.krw0_logedit.text()), 0)
            assert krw0 >= 0.
            krw0 = krw0 * 10. ** krw0_log
            if not self.mod_krw:
                params['krw0'][1] = krw0
            else:
                self.parent().flood.krw = krw0
            # self.update_data_lines()

        except Exception as e:
            msg = QMessageBox(parent=self, text=str(e))
            msg.exec_()
            if not self.mod_krw:
                krw0 = params['krw0'][1]
            else:
                krw0 = self.parent().flood.krw

        self.silent_set_params(params['C'][1], krw0)
        self.update_data_lines()
        self.set_perm_text()
        self.parent().flood.permeability[0] = krw0 * params['perm'][1]

    def silent_set_params(self, c: float, krw0: float):

        self.c_edit.blockSignals(True)
        self.c_edit.setText(str(c))
        self.c_edit.blockSignals(False)

        self.krw0_edit.blockSignals(True)
        self.krw0_logedit.blockSignals(True)
        log_level = np.floor(np.log10(krw0))
        self.krw0_edit.setText(str(np.round(krw0 * 10. ** -log_level, 1)))
        self.krw0_logedit.setText(str(np.round(log_level)))
        self.krw0_edit.blockSignals(False)
        self.krw0_logedit.blockSignals(False)

    def update_data_lines(self):

        for i in range(len(self.data_lines)):
            line = self.data_lines.pop()
            try:
                self.chart.removeSeries(line)
            except Exception as e:
                print(e)
            del line

        fv = self.parent()
        p_flow_rates, plateaus = fv.get_rates_and_plateaus(self.i)
        offset = fv.flood.used_offsets[self.i]
        plateaus = np.array(plateaus) - offset
        print(p_flow_rates, plateaus)
        my_flood = fv.flood
        # my_fluid = my_flood.fluid

        my_core = my_flood.core
        # us = np.divide(np.divide(p_flow_rates, 60.), my_core.my_area())

        my_experiment = fv.linked_widget.experiment
        shear = my_flood.get_plateau_shear_rates(my_experiment, [])
        # saturation = my_experiment.get_flood_saturation(my_flood)
        # c = my_experiment.petro_parameters['C'][1]
        k = my_experiment.petro_parameters['perm'][1]
        if not self.mod_krw:
            krw0 = my_experiment.petro_parameters['krw0'][1]
        else:
            krw0 = self.parent().flood.krw
        # phi = my_experiment.petro_parameters['phi'][1]
        # n = my_fluid.specific_fluid.polymer_solution.model['n']
        # shear = c * np.divide(us, np.sqrt(0.5 * (krw0 * k * 0.0000009869 / 1000.) * phi * saturation))
        # shear = ((3. * n + 1.) / (4. * n)) ** (n / (n - 1.)) * shear
        print('Shear rates:', shear)
        mu = np.multiply(krw0 * k * my_core.my_area() / (245. * my_core.my_length()), np.divide(plateaus, p_flow_rates))
        print('Polymer viscosity:', mu)

        for s, m in zip(shear, mu):
            data_line = QScatterSeries()
            data_line.setColor(self.color)
            data_line.setMarkerShape(QScatterSeries.MarkerShapeCircle)
            data_line.setMarkerSize(10.)
            data_line.setUseOpenGL(False)
            self.chart.addSeries(data_line)
            data_line.attachAxis(self.x_axis)
            data_line.attachAxis(self.y_axis)
            point = series_to_polyline([s], [m])
            data_line.append(point)
            self.data_lines.append(data_line)

    def mousePressEvent(self, a0: QtGui.QMouseEvent):
        print('PolymerPermView mouse press.', self, a0)

    def keyPressEvent(self, a0: QtGui.QKeyEvent):

        pps = self.parent().linked_widget.experiment.petro_parameters
        # c = np.round(float(self.c_edit.text()), 1)
        c = pps['C'][1]
        if not self.mod_krw:
            krw0 = pps['krw0'][1]
        else:
            krw0 = self.parent().flood.krw

        if a0.key() == 52:
            if c - 0.1 >= 1.:
                c -= 0.1
                self.silent_set_params(c, krw0)
                self.update_params()
        elif a0.key() == 54:
            c += 0.1
            self.silent_set_params(c, krw0)
            self.update_params()
        elif a0.key() == 56:
            log_level = np.floor(np.log10(krw0))
            pre_k = krw0 * 10. ** -log_level + 0.1
            self.silent_set_params(c, pre_k * 10. ** log_level)
            self.update_params()
        elif a0.key() == 50:
            log_level = np.floor(np.log10(krw0))
            pre_k = krw0 * 10. ** -log_level - 0.1
            self.silent_set_params(c, pre_k * 10. ** log_level)
            self.update_params()
        else:
            print(a0.key(), Qt.Key_Up)

    def closeEvent(self, a0: QtGui.QCloseEvent):

        try:
            pps = self.parent().linked_widget.experiment.petro_parameters
            if pps['C'][2] is None:
                pps['C'][1] = pps['C'][0]
            if pps['krw0'][2] is None:
                pps['krw0'][1] = pps['krw0'][0]

            self.parent().perm_view = None

        except Exception as e:
            print(e)

    def set_c(self):
        pps = self.parent().linked_widget.experiment.petro_parameters
        pps['C'][2] = self.parent().flood

    def set_krw0(self):
        if not self.mod_krw:
            pps = self.parent().linked_widget.experiment.petro_parameters
            pps['krw0'][2] = self.parent().flood


class OffsetTool(QDialog):

    def __init__(self, parent):
        super(OffsetTool, self).__init__(parent=parent)

        # self.setFixedSize(100, 250)

        try:
            self.current_offsets = parent.flood.offsets
            self.current_manual_offsets = parent.flood.manual_offsets
            self.current_used_offsets = parent.flood.used_offsets
            self.current_offset_source = parent.flood.offset_source_flood
        except Exception as e:
            msg = QMessageBox(parent=self, text=str(e))
            msg.exec_()
            return

        layout = QGridLayout()
        self.setLayout(layout)
        self.edits = []

        for i in range(5):
            if i == 0:
                txt = 'W'
            else:
                txt = str(i)

            layout.addWidget(QLabel(parent=self, text=txt + ' [psi]:'), i, 0)
            self.edits.append(QLineEdit(parent=self))
            self.edits[i].setText(str(parent.flood.manual_offsets[i]))
            self.edits[i].editingFinished.connect(self.check_values)
            layout.addWidget(self.edits[i], i, 1)

        try:
            self.reset_button = QPushButton(parent=self, text='Reset')
            self.reset_button.clicked.connect(lambda: self.update_view_offsets(True))
            layout.addWidget(self.reset_button, 6, 0)
        except Exception as e:
            msg = QMessageBox(parent=self, text=str(e))
            msg.exec_()
            return

        self.show()

    def check_values(self):
        print('calling check values')

        is_error = False

        for i, edit in enumerate(self.edits):
            try:
                float(edit.text())

            except ValueError:
                is_error = True
                edit.setText(str(self.parent().flood.manual_offsets[i]))

        if not is_error:
            self.update_view_offsets()

    def update_view_offsets(self, reset_offsets: bool=False):
        print('calling update view offsets')

        fv = self.parent()
        fl = fv.flood

        if not reset_offsets:

            for i, edit in enumerate(self.edits):
                fl.manual_offsets[i] = float(edit.text())

            fl.used_offsets = fl.manual_offsets
            fl.offset_source_flood = fl

        else:

            for i, edit in enumerate(self.edits):
                edit.setText(str(self.current_manual_offsets[i]))

            fl.manual_offsets = self.current_manual_offsets
            fl.offsets = self.current_offsets
            fl.used_offsets = self.current_used_offsets
            fl.offset_source_flood = self.current_offset_source

        for i in range(5):
            td = fl.time_data.copy()
            pd = fl.p_data[i].copy()
            fl.set_curve_data(td, pd, i)

        if fv.mode == 'Zoom':
            fl.set_x_data_view_time()
            fv.zoom_mode()
        elif fv.mode == 'Plateau':
            fl.set_x_data_view_PV()
            fv.plateau_mode()
        else:
            fl.set_x_data_view_time()
            fv.flow_mode()

    def closeEvent(self, a0: QtGui.QCloseEvent):

        for edit in self.edits:
            edit.blockSignals(True)

        fv = self.parent()

        try:
            fv.offsets_view = None
        except Exception as e:
            print(e)


if __name__ == "__main__":

    app = QApplication(sys.argv)
    f = Fluid()
    cvf = CoreViewForm()
    c = cvf.core_view.core
    form = CoreFloodForm(f, c)
    form.move(0, 0)
    form.show()
    cvf.current_flood_form = form
    sys.exit(app.exec_())
