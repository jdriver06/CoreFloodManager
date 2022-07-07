
from PyQt5.Qt import QPolygonF, QFileDialog, pyqtSignal, Qt
from PyQt5.QtWidgets import QMessageBox, QListWidget, QListWidgetItem, QVBoxLayout, QHBoxLayout, \
    QWidget, QGroupBox, QLabel, QPushButton, QTableWidget, QTableWidgetItem
from PyQt5.QtGui import QMouseEvent, QKeyEvent
# from PyQt5.QtSql import QSqlQuery
import numpy as np
from enum import Enum, auto


__version__ = '0.1.0'


def list_element_wise_compare(lst: list, element) -> list:

    bool_list = []
    for item in lst:
        try:
            cmp = item == element
        except Exception as e:
            cmp = False
            print(e)
        bool_list.append(cmp)

    return bool_list


def series_to_polyline(x, y):

    size = len(x)
    polyline = QPolygonF(size)
    pointer = polyline.data()
    d_type, t_info = np.float, np.finfo
    pointer.setsize(2 * polyline.size() * t_info(d_type).dtype.itemsize)
    memory = np.frombuffer(pointer, d_type)
    memory[:(size - 1) * 2 + 1:2] = x
    memory[1:(size - 1) * 2 + 2:2] = y

    return polyline


def find_order_above(num, start=0., base=10.):

    order = start
    found = False
    while not found:
        mod = divmod(num, base ** order)
        if mod[0] == 0.:
            found = True
        else:
            order += 1.

    return order


def find_order_below(num, start=0., base=10.):

    order = start
    found = False
    while not found:
        mod = divmod(base ** order, num)
        if mod[0] == 0.:
            found = True
        else:
            order -= 1.

    return order


def find_nearest(x_vector: np.ndarray, y_vector: np.ndarray, x, y, mx: float=1., my: float=1.):

    if not x_vector.size or not y_vector.size:
        return np.NaN

    d = np.power((x_vector - x) * mx, 2.) + np.power((y_vector - y) * my, 2.)
    index = np.where(d == d.min())

    return index[0][0]


def file_open_dlg(direct, filter, save=False):

    options = QFileDialog.Options()
    options |= QFileDialog.DontUseNativeDialog

    if save:
        file_name, _ = QFileDialog.getSaveFileName(None, "QFileDialog.getSaveFileName()",
                                                   direct, filter, options=options)
    else:
        file_name, _ = QFileDialog.getOpenFileName(None, "QFileDialog.getOpenFileName()",
                                                   direct, filter, options=options)

    return file_name


# class JQSqlQuery(QSqlQuery):
#
#     def __init__(self, **kwargs):
#         super(JQSqlQuery, self).__init__(**kwargs)
#
#     def select(self, table_name: str, columns: list, execute: bool):
#
#         query_string = 'SELECT '
#
#         for i, column in enumerate(columns):
#             query_string += column
#             if i < len(columns) - 1:
#                 query_string += ', '
#             else:
#                 query_string += ' '
#
#         query_string += 'FROM ' + table_name
#
#         if execute:
#             self.exec_(query_string)
#         else:
#             self.prepare(query_string)


class SignalList:

    def __init__(self, parent, i: int):

        self.parent = parent
        self.index = i
        self.ref_lists = []
        self.referencers = []
        self.ref_parents = []
        self.ref_indecies = []

        self.current_row = -1

        self.item_names = []
        self.objects = []

    def set_ref_lists(self, ref_lists: list):

        self.ref_lists = ref_lists
        self.ref_parents = []
        self.ref_indecies = []

        for c_list in ref_lists:
            c_list.add_referencer(self)
            self.ref_parents.append(c_list.parent)
            index = -1
            for i, cc_list in enumerate(self.ref_parents[-1].signal_lists):
                if cc_list == c_list:
                    index = i
                    break
            self.ref_indecies.append(index)

    def check_name_before_adding(self, name):

        if name in self.item_names:
            return True

        return False

    def check_item_before_removing(self, obj):

        for c_obj in self.objects:
            if hasattr(c_obj, 'ref_objects'):
                for r_obj in c_obj.ref_objects:
                    if obj == r_obj:
                        return True

        return False

    def add_referencer(self, referencer):

        if not isinstance(referencer, SignalList):
            raise TypeError('SignalList can only be referenced by other SignalList objects.')

        if referencer not in self.referencers:
            self.referencers.append(referencer)

    def find_item_by_name(self, name: str):

        for obj in self.objects:
            if obj.name == name:
                return obj

        return None


class SignalListWidget(QListWidget):

    selected = pyqtSignal(int)

    def __init__(self, parent, sl: SignalList):
        super(SignalListWidget, self).__init__(parent=parent)
        font = self.font()
        font.setPointSize(10)
        self.setFont(font)
        self.signal_list = sl
        for name in sl.item_names:
            self.addItem(QListWidgetItem(name), False)
        self.itemClicked.connect(self.update_current_row_in_list)

    def set_ref_lists(self, ref_lists: list):
        self.signal_list.ref_lists = ref_lists
        self.signal_list.ref_parents = []
        self.signal_list.ref_indecies = []
        for c_list in ref_lists:
            c_list.referencers.append(self)
            self.signal_list.ref_parents.append(c_list.parent())
            index = -1
            for i, cc_list in enumerate(self.signal_list.ref_parents[-1].lists):
                if cc_list == c_list:
                    index = i
                    break
            self.signal_list.ref_indecies.append(index)

    def mousePressEvent(self, e: QMouseEvent):
        super(SignalListWidget, self).mousePressEvent(e)
        self.selected.emit(self.signal_list.index)

    def takeItem(self, row: int):

        obj = self.signal_list.objects[row]
        for referencer in self.signal_list.referencers:
            if referencer.check_item_before_removing(obj):
                return None

        item = super(SignalListWidget, self).item(row)
        super(SignalListWidget, self).takeItem(row)
        self.signal_list.item_names.pop(row)
        if hasattr(obj, '__del__'):
            obj.__del__()

        self.signal_list.objects.pop(row)

        return item

    def addItem(self, aitem: QListWidgetItem, add_to_sl: bool=True):

        if self.signal_list.check_name_before_adding(str(aitem)) and add_to_sl:
            return False

        super(SignalListWidget, self).addItem(aitem)

        if add_to_sl:
            self.signal_list.item_names.append(str(aitem))

        return True

    def update_current_row_in_list(self):
        self.signal_list.current_row = self.currentRow()


class SignalListManager:

    def __init__(self):
        self.signal_lists = []

    def add_signal_list(self, n: int=1):
        for i in range(n):
            count = len(self.signal_lists)
            self.signal_lists.append(SignalList(self, count))


class SignalListManagerWidget(QWidget):

    def __init__(self, parent, pml: SignalListManager, names: list, cls: list, *args):
        super(SignalListManagerWidget, self).__init__(parent=parent)

        self.setWindowTitle(' ')

        self.cls = cls
        self.args = args
        self.pm_list = pml
        # self.objects = []

        # if len(cls[1]) > 0:
        #     for i in range(len(cls[1])):
        #         self.objects.append([])

        self.names = names
        self.selected_list = 0

        lw_layout = QVBoxLayout()
        self.setLayout(lw_layout)

        lw_group_box = QGroupBox()
        lgb_layout = QHBoxLayout()
        lw_group_box.setLayout(lgb_layout)

        self.lists = []
        list_layout = QHBoxLayout(self)
        for i in range(len(names)):
            sub_layout = QVBoxLayout(self)
            text = QLabel(parent=self, text=names[i] + ':')
            font = text.font()
            font.setPointSize(10)
            # font.setBold(True)
            text.setFont(font)
            sub_layout.addWidget(text)
            list_layout.addLayout(sub_layout)
            self.lists.append(SignalListWidget(self, self.pm_list.signal_lists[i]))
            self.lists[i].selected.connect(self.select_list)
            self.lists[i].doubleClicked.connect(self.view_item)
            sub_layout.addWidget(self.lists[i])

        self.add_button = QPushButton(parent=self, text='Add')
        self.add_button.clicked.connect(self.add_item)
        self.remove_button = QPushButton(parent=self, text='Remove')
        self.remove_button.clicked.connect(self.remove_item)

        font = self.add_button.font()
        font.setPointSize(10)
        self.add_button.setFont(font)
        self.remove_button.setFont(font)

        lgb_layout.addWidget(self.add_button)
        lgb_layout.addWidget(self.remove_button)

        lw_layout.addLayout(list_layout)
        lw_layout.addWidget(lw_group_box)

        self.lists[0].selected.emit(0)

    def populate_from_pm_list(self):

        for signal_list_widget, signal_list in zip(self.lists, self.pm_list.signal_lists):
            signal_list_widget.signal_list = signal_list
            for name in signal_list.item_names:
                item = QListWidgetItem(name)
                signal_list_widget.addItem(item, False)

    def add_item(self):
        cls = self.cls[0][self.selected_list]
        args = self.args
        print(cls, args)
        try:
            cls(self, *args)
        except Exception as e:
            print('Error with add item: ', e)

    def remove_item(self):

        i = self.lists[self.selected_list].currentRow()
        item = self.lists[self.selected_list].takeItem(i)

        if item is None:
            msg = QMessageBox(parent=self.parent(), text='Cannot delete item in use.')
            msg.exec_()
            return

        del item

    def view_item(self):
        i = self.lists[self.selected_list].currentRow()
        name = self.lists[self.selected_list].item(i)
        name = name.text()
        obj = self.pm_list.signal_lists[self.selected_list].objects[i]

        try:
            print('view_item', obj, self.parent())
            obj.view_class(obj, self.parent())
        except Exception as e:
            print('Error with view item:', e)
            try:
                obj.view_class(self, obj, name, False)
            except Exception as e2:
                print(e2)
                try:
                    obj.view_class(obj, self.parent(), [])
                except Exception as e3:
                    print(e3)
        finally:
            return

    def select_list(self, i: int):

        self.selected_list = i
        for j in range(len(self.lists)):
            if j == i:
                self.lists[j].setStyleSheet('border: 1px solid black')
            else:
                self.lists[j].setStyleSheet('border: 1px solid lightgray')

    def submit_name(self, name: str, *args):

        s_list = self.pm_list.signal_lists[self.selected_list]
        s_view = self.lists[self.selected_list]

        if self.cls[1]:
            cls = self.cls[1][self.selected_list]
            ref_objs = []
            if cls != str:
                if s_list.ref_lists:
                    for ref_list in s_list.ref_lists:
                        ref_obj = ref_list.current_row
                        if ref_obj == -1:
                            return
                        ref_obj = ref_list.objects[ref_obj]
                        ref_objs.append(ref_obj)
                    obj = cls(name, *ref_objs, *args)
                    obj.ref_objects = ref_objs
                else:
                    obj = cls(name, *args)
                obj.project_manager_list = self
            else:
                obj = None

            if s_view.addItem(name):
                s_list.objects.append(obj)

    def submit_name_ref_objects(self, name: str, ref_objs: list, *args):

        print('submit_name_ref_objects called.', self.selected_list, self.pm_list)

        s_list = self.pm_list.signal_lists[self.selected_list]
        s_view = self.lists[self.selected_list]

        print(self.cls[1][self.selected_list])

        if self.cls[1]:
            cls = self.cls[1][self.selected_list]
            if cls != str:
                if s_list.ref_lists:
                    try:
                        # why is this done this way? the exception case is more natural.
                        print('first try:')
                        print(cls, name, ref_objs, args)
                        obj = cls(name, *ref_objs, *args)
                        obj.ref_objects = ref_objs
                    except Exception as e:
                        print('error with submit_name_ref_objects, using contingency')
                        obj = cls(name, ref_objs, *args)
                        print('object created:', obj)
                else:
                    obj = cls(name, *args)

                obj.project_manager_list = self

            else:
                obj = None

            try:
                if s_view.addItem(name):
                    s_list.objects.append(obj)
            except Exception as e:
                print(e)

            print('submit name ref objects finished.')

    def clear_items(self):

        for i in range(len(self.names)):
            self.pm_list.signal_lists[i].objects = []
            count = self.lists[i].count()
            for j in range(count):
                item = self.lists[i].takeItem(self.lists[i].item(count-1-j))
                del item


class SignalTable(QTableWidget):

    rc_changed = pyqtSignal(int)
    data_changed = pyqtSignal(int, int, float)
    rows_changed = pyqtSignal(int)

    class ThresholdMode(Enum):

        GREATER_THAN = auto()
        GREATER_THAN_OR_EQUAL_TO = auto()
        LESS_THAN = auto()
        LESS_THAN_OR_EQUAL_TO = auto()

    def __init__(self, parent, data: np.array):
        super(SignalTable, self).__init__(parent=parent)

        self.data = data
        self.setRowCount(data.shape[1] + 1)
        self.setColumnCount(data.shape[0])
        self.threshold_modes = []
        self.thresholds = []

        self.initialize_thresholds_and_modes()
        self.initialize_table()

        self.itemChanged.connect(self.check_item)
        self.rc_changed.connect(self.update_rows)

    def initialize_thresholds_and_modes(self):

        t_modes = []
        ts = []
        for c in range(self.columnCount()):
            t_modes.append(self.ThresholdMode.GREATER_THAN)
            ts.append(np.nan)

        self.threshold_modes = t_modes
        self.thresholds = ts

    def sync_item(self, r: int, c: int):

        item = self.item(r, c)
        if item is None:
            item = QTableWidgetItem()
            self.setItem(r, c, item)

        self.blockSignals(True)
        if not np.isnan(self.data[c, r]):
            item.setText(str(self.data[c, r]))
        else:
            item.setText('')
        self.blockSignals(False)

    def initialize_table(self):

        for r in range(self.data.shape[1]):
            for c in range(self.data.shape[0]):
                self.sync_item(r, c)

    def check_val(self, val: float, c: int) -> bool:

        if self.threshold_modes[c] == self.ThresholdMode.GREATER_THAN:
            if val > self.thresholds[c]:
                return True
            return False

        elif self.threshold_modes[c] == self.ThresholdMode.GREATER_THAN_OR_EQUAL_TO:
            if val >= self.thresholds[c]:
                return True
            return False

        elif self.threshold_modes[c] == self.ThresholdMode.LESS_THAN:
            if val < self.thresholds[c]:
                return True
            return False

        else:
            if val <= self.thresholds[c]:
                return True
            return False

    def check_item(self, item: QTableWidgetItem):

        txt = item.text()
        r = item.row()
        c = item.column()

        try:

            val = float(txt)
            if not self.check_val(val=val, c=c):
                raise ValueError('Value outside threshold.')

            if r == self.rowCount() - 1:
                self.update_rows(self.rowCount() + 1)
            else:
                self.data[c, r] = val

            self.data_changed.emit(r, c, val)

        except ValueError as e:

            print('ValueError:', e)
            if r < self.rowCount() - 1:
                self.sync_item(r, c)
            else:
                self.takeItem(r, c)
                del item

    def keyPressEvent(self, e: QKeyEvent):
        super(SignalTable, self).keyPressEvent(e)

        if e.key() == Qt.Key_Delete:

            r = self.currentRow()
            c = self.currentColumn()

            if r == self.rowCount() - 1:
                self.update_rows(self.rowCount() - 1)

            else:
                self.data[c, r] = np.nan
                self.data_changed.emit(r, c, np.nan)
                self.sync_item(r, c)

    def update_rows(self, nrc: int):

        rc = self.rowCount()

        if nrc < rc:

            self.data = self.data[:, :nrc - 1]

            for c in range(self.columnCount()):
                item = self.takeItem(nrc - 1, c)
                del item

            self.rows_changed.emit(nrc)

        else:

            new_data = np.array(np.zeros((self.columnCount(), nrc - 1)))
            new_data[:, :rc - 1] = self.data[:, :]
            self.rows_changed.emit(nrc)

            for r in range(rc, nrc):
                for c in range(self.columnCount()):

                    item = self.item(r - 1, c)

                    if item is not None:
                        val = float(item.text())
                        new_data[c, r - 1] = val
                        self.data_changed.emit(r, c, val)
                    else:
                        new_data[c, r - 1] = np.nan

            self.data = new_data

        self.setRowCount(nrc)
