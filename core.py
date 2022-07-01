
from __future__ import division
from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot, Qt, QEvent
from PyQt5.QtWidgets import QDialog, QLineEdit, QLabel, QComboBox, QWidget, QVBoxLayout, \
    QMenu, QAction, QMainWindow, QFileDialog, QTextEdit, QListWidget, QListWidgetItem, QMessageBox, QCheckBox
from PyQt5.Qt import QRect, QApplication, QPainter, QIcon, QBrush, QPoint, QPen, QColor, QPushButton, QPointF
from PyQt5.QtGui import QPaintEvent, QMouseEvent, QFocusEvent
import numpy as np
import copy
import pickle as pkl


__version__ = '0.1.0'


class Core:

    def __init__(self, name: str,  c_sec, diameters=(), lengths=(), masses=(), lithology='sandstone', pm=None):
        self.name = name
        self.project_manager_list = pm
        self.view_class = CoreViewForm
        self.core_sections = []
        self.cs_rects = []
        self.ref_objects = []
        self.cut_refs = [[]]
        self.reverse_loaded = [False]
        self.add_existing_section(c_sec)

        for d, l, m in zip(diameters, lengths, masses):
            cs = CoreSection('dummy_name', d, l, m, lithology)
            cs.project_manager = pm
            self.core_sections.append(cs)
            r = QRect(50, 50, 50, 50)
            r.setWidth(int(20*d))
            r.setHeight(int(20*l))

            self.cs_rects.append(r)
            self.cut_refs[0].append(0)

    def add_existing_section(self, cs):

        if cs in self.core_sections:
            return False

        if len(cs.cuts) > 1:
            core_list = self.project_manager_list.pm_list.signal_lists[0]
            i = core_list.objects.index(self)
            if i + 1 < len(core_list.objects):
                for j in range(i + 1, len(core_list.objects)):
                    c = core_list.objects[j]
                    if cs in c.core_sections:
                        return False

        d = cs.diameter
        cs.length = cs.cuts[-1]
        l = cs.length

        r = QRect(50, 50, 50, 50)
        r.setWidth(int(20 * d))
        r.setHeight(int(20 * l))
        # t_h = 0
        # for rect in self.cs_rects:
        #     t_h += rect.height()
        # r.moveTo(QPoint(50, 640 - r.height() - t_h))

        self.core_sections.append(cs)
        if cs not in self.ref_objects:
            self.ref_objects.append(cs)
        self.cs_rects.append(r)
        self.set_rect_positions()
        self.add_to_cut_refs(len(cs.cuts) - 1)

        return True

    def add_section(self, d, l, m, lithology=''):

        if lithology == '':
            lithology = self.core_sections[0].lithology

        cs = CoreSection('dummy_name', d, l, m, lithology)

        r = QRect(50, 50, 50, 50)
        r.setWidth(int(20*d))
        r.setHeight(int(20*l))
        # t_h = 0
        # for rect in self.cs_rects:
        #     t_h += rect.height()
        # r.moveTo(QPoint(50, 640 - r.height() - t_h))
        self.core_sections.append(cs)
        self.cs_rects.append(r)
        self.set_rect_positions()
        self.add_to_cut_refs()

    def set_rect_positions(self, cut_num: int=0):

        # temp = self.cs_rects.copy()
        if self.reverse_loaded[cut_num]:
            cs_rects = reversed(self.cs_rects)
        else:
            cs_rects = self.cs_rects

        t_h = 0
        for rect in cs_rects:
            rect.moveTo(QPoint(50, 640 - rect.height() - t_h))
            t_h += rect.height()

    def add_to_cut_refs(self, cut_num: int=0):

        for cut_ref in self.cut_refs:
            cut_ref.append(cut_num)

    def remove_section(self, cut_num: int=0):

        if len(self.cs_rects) < 1:
            return

        self.cs_rects.pop()
        # self.reverse_loaded.pop()
        for cut_ref in self.cut_refs:
            cut_ref.pop()

        cs = self.core_sections.pop()

        if cs in self.ref_objects:
            self.ref_objects.remove(cs)

        self.set_rect_positions(cut_num)

    def set_section_attributes(self, i: int, d: float, l: float, m: float, lithology: str,
                               phi: float, notes: str, cut_num: int=-1):

        cs = self.core_sections[i]
        r = self.cs_rects[i]
        cs.diameter = d
        cs.length = l
        # cs.cuts[self.cut_refs[cut_num][i]] = l
        cs.mass = m
        cs.lithology = lithology
        cs.est_porosity = phi
        cs.notes = notes
        r.setWidth(int(20*d))
        r.setHeight(int(20*l))

    def create_cut(self, lengths: list, rev: bool=False):

        if len(lengths) != len(self.core_sections):
            return

        cr = self.cut_refs[-1]
        ncr = []
        cores = self.project_manager_list.pm_list.signal_lists[0].objects
        ind = cores.index(self)
        if ind + 1 < len(cores):
            cores = cores[ind + 1:]

        i = -1
        for length, section in zip(lengths, self.core_sections):
            i += 1

            if length < section.cuts[cr[i]]:
                in_use = False
                for c in cores:
                    if section in c.core_sections:
                        in_use = True
                if not in_use:
                    section.cut(length)
                    ncr.append(len(section.cuts) - 1)
                else:
                    ncr.append(cr[i])
                    msg = QMessageBox(parent=self.project_manager_list.parent(),
                                      text='Section #' + str(i + 1) +
                                           ' cannot be cut because it is in use by a subsequent core.')
                    msg.exec_()
            else:
                ncr.append(cr[i])

        self.cut_refs.append(ncr)
        self.reverse_loaded.append(rev)

    def my_estimated_pore_volume(self):
        pv = 0.0
        for section in self.core_sections:
            pv += section.my_estimated_pore_volume()

        return pv

    def my_length(self):

        L = 0.

        for section in self.core_sections:
            L += section.length

        return L

    def my_area(self):

        A = []
        lens = []

        for section in self.core_sections:
            A.append(section.my_area())
            lens.append(section.length)

        A = np.sum(np.multiply(A, lens)) / np.sum(lens)
        # sA = sum(A)
        # l = float(len(A))
        #
        # A = sA / l

        return A

    def my_bulk_volume(self):

        bv = 0.

        for section in self.core_sections:
            bv += section.my_area() * section.length

        return bv

    def load_from_pickle(self, pickle):

        l = len(pickle.diameters)
        l_cs = len(self.core_sections)

        self.set_section_attributes(0, pickle.diameters[0], pickle.lengths[0], pickle.masses[0],
                                    pickle.lithologies[0], pickle.phis[0], pickle.notes[0])

        if l > l_cs:
            for i in range(0, l - l_cs):
                self.add_section(1, 1, 1)
        elif l < l_cs:
            for i in range(0, l_cs - l):
                self.remove_section()

        i = 0
        for d, l, m, lithology, phi, note in zip(pickle.diameters, pickle.lengths,
                                               pickle.masses, pickle.lithologies, pickle.phis, pickle.notes):

            self.set_section_attributes(i, d, l, m, lithology, phi, note)
            i += 1

    def construct_pickle(self):

        d = []
        l = []
        m = []
        liths = []
        phis = []
        notes = []

        for cs in self.core_sections:
            d.append(cs.diameter)
            l.append(cs.length)
            m.append(cs.mass)
            liths.append(cs.lithology)
            phis.append(cs.est_porosity)
            notes.append(cs.notes)

        return CoreSaveStructure(d, l, m, liths, phis, notes)

    def __del__(self):

        # print('deleting core.', self)

        cores = self.project_manager_list.pm_list.signal_lists[0].objects
        if self not in cores:
            return

        for i, section in enumerate(self.core_sections):
            cr = []
            for cut_ref in self.cut_refs:
                cr.append(cut_ref[i])

            keep = []
            for cut in section.cuts:
                keep.append(False)
            keep[0] = True

            for c in cores:
                if c == self:
                    pass
                else:
                    if section in c.core_sections:
                        j = c.core_sections.index(section)
                        for cut_ref in c.cut_refs:
                            keep[cut_ref[j]] = True

            ncr = [j for j, k in enumerate(keep) if k]

            for c in cores:
                if c == self:
                    pass
                else:
                    if section in c.core_sections:
                        j = c.core_sections.index(section)
                        # print(j, c.cut_refs, ncr, keep)
                        for cut_ref in c.cut_refs:
                            cut_ref[j] = ncr.index(cut_ref[j])

            new_cuts = [cut for k, cut in zip(keep, section.cuts) if k]
            # print(new_cuts)
            section.cuts = new_cuts


class CoreSaveStructure:

    def __init__(self, diameters: list, lengths: list, masses: list, lithologies: list, phis: list, notes: list):
        self.diameters = diameters
        self.lengths = lengths
        self.masses = masses
        self.lithologies = lithologies
        self.phis = phis
        self.notes = notes


class CoreSection:

    # a dictionary of matrix densities
    densities = {"sandstone": 2.65, "limestone": 2.71, "dolomite": 2.87}
    # signal to recalculate derived properties
    # measurementChanged = pyqtSignal()

    def __init__(self, name: str, diameter, length, mass, lithology, phi=0.22, notes='', pm=None):
        # super(CoreSection, self).__init__()
        self.name = name
        self.project_manager_list = pm
        self.diameter = diameter  # section diameter in cm
        self.length = length  # section length in cm
        self.cuts = [length]
        self.mass = mass  # section weight in grams
        self.est_porosity = phi
        self.notes = notes
        if lithology in CoreSection.densities:
            self.lithology = lithology
        # self.measurementChanged.connect(self.mc)
        # self.measurementChanged.emit()

    @staticmethod
    def calc_area(diameter):
        return 0.25*np.pi*diameter**2.0

    def my_area(self):
        return CoreSection.calc_area(self.diameter)

    def my_bulk_volume(self, cut_num: int=0):
        return self.my_area()*self.cuts[cut_num]

    def my_estimated_porosity(self):
        return 1.0-(self.mass/self.my_bulk_volume())/CoreSection.densities[self.lithology]

    def my_estimated_pore_volume(self):
        phi = self.my_estimated_porosity()
        BV = self.my_bulk_volume()
        return phi * BV

    def my_estimated_fluid_density(self):
        rho_matrix = CoreSection.densities[self.lithology]
        return (self.mass/self.my_bulk_volume() - rho_matrix * (1. - self.est_porosity)) / self.est_porosity

    def cut(self, l: float):

        if l > self.cuts[-1] or l < 0:
            return

        self.cuts.append(l)


class CoreHolder:

    def __init__(self, taps=(7.62, 15.24, 22.86), confining_pressure=1000.0):
        self.taps = taps  # tap positions in cm
        self.confining_pressure = confining_pressure  # confining pressure in psi
        self.rect = QRect(20, 20, 135, 630)


class CoreViewForm(QMainWindow):

    def __init__(self, c: Core=None, parent=None):

        super(CoreViewForm, self).__init__(parent=parent)
        self.core_view = CoreView(self, c)
        self.current_flood_form = None
        self.setCentralWidget(self.core_view)

        if c is not None:
            self.setWindowTitle('Core Viewer: ' + c.name)
        else:
            self.setWindowTitle('Core Viewer: New_Core')
        icon = QIcon('Coreholder Cropped.jpg')
        self.setWindowIcon(icon)

        self.file_menu = QMenu('&File', self)
        self.file_menu.addAction('&New', self.file_open)
        self.file_menu.addAction('&Open', self.file_open)
        self.file_menu.addAction('&Save', self.file_save)
        self.menuBar().addMenu(self.file_menu)

        self.setFixedSize(175, 690)

        print('{} Length: {:.2f} cm, Area: {:.2f} cm2, BV: {:.1f} cm3'.format(c.name, c.my_length(), c.my_area(),
                                                                              c.my_bulk_volume()))

        self.show()

    def size_to_widget(self):
        self.setFixedSize(self.core_view.width(), self.core_view.height() + 20)

    def file_open(self):

        dlg = QFileDialog()
        options = dlg.options()
        options |= dlg.DontUseNativeDialog
        dlg.setOptions(options)
        file_name, _ = dlg.getOpenFileName(self, "Select core", "",
                                           "Pickle Files (*.core)")
        if file_name:
            self.load_pickle(file_name)
            ff = self.current_flood_form
            if ff is not None:
                ff.flood.core = self.core_view.core
                ff.update_plateau_line_x()
                if ff.x_view == 'PV':
                    ff.plateau_mode()

    def load_pickle(self, file_name: str):

        pickle_in = open(file_name, 'rb')
        pickle = pkl.load(pickle_in)
        pickle_in.close()
        self.core_view.core.load_from_pickle(pickle)
        self.core_view.set_section_positions()
        i = file_name.find('.')
        j = file_name.rfind('/')
        core_name = file_name[j + 1:i]
        self.setWindowTitle('Core Viewer: ' + core_name)

    def file_save(self):

        CoreSaveDlg(self)


class CoreSaveDlg(QDialog):

    def __init__(self, parent: CoreViewForm):

        super(CoreSaveDlg, self).__init__(parent)

        self.setWindowTitle('Save Core')

        lyt = QVBoxLayout()
        self.parent = parent
        self.name_label = QLabel('Core Name:')
        self.name_edit = QLineEdit('')
        self.submit_button = QPushButton('Save')
        self.submit_button.clicked.connect(self.save_core)

        lyt.addWidget(self.name_label)
        lyt.addWidget(self.name_edit)
        lyt.addWidget(self.submit_button)
        self.setLayout(lyt)

        self.show()

    def save_core(self):

        pickle = self.parent.core_view.core.construct_pickle()
        name = self.name_edit.text()
        pickle_out = open(name + '.pickle', 'wb')
        pkl.dump(pickle, pickle_out)
        pickle_out.close()
        self.close()


class PlugSelector(QDialog):

    def __init__(self, parent):
        super(PlugSelector, self).__init__(parent=parent)

        layout = QVBoxLayout(self)
        self.setLayout(layout)

        self.label = QLabel(parent=self, text='Plug List:')
        self.list = QListWidget(self)
        layout.addWidget(self.label)
        layout.addWidget(self.list)

        self.list.itemDoubleClicked.connect(self.add_item_to_core)

    def add_item(self, item: QListWidgetItem):

        self.list.addItem(item)

    def add_item_to_core(self):

        i = self.list.currentIndex().row()
        plw = self.parent().core.project_manager_list.parent().plug_list_widget
        cs = plw.pm_list.signal_lists[0].objects[i]
        is_added = self.parent().core.add_existing_section(cs)
        if not is_added:
            msg = QMessageBox(parent=self, text='Plug cannot be added to core.')
            msg.exec()
        else:
            self.parent().core.add_existing_section(cs)
            # self.parent().core.ref_objects.append(cs)
            self.parent().core.set_rect_positions(self.parent().cut_selector.currentIndex())


class CoreView(QWidget):

    def __init__(self, cvf: CoreViewForm, c: Core):

        super(CoreView, self).__init__()

        self.form = cvf
        self.setFixedSize(175, 670)

        if c is not None:
            self.core = c
        else:
            self.core = Core([], [], [], 'sandstone')

        self.core_holder = CoreHolder()

        self.core_entry_form = CoreEntryForm(self)

        self.cut_text = QLabel(parent=self, text='Cut Number:')
        self.cut_text.setFixedSize(80, 25)
        self.cut_text.move(QPoint(180, 5))
        self.cut_selector = QComboBox(parent=self)
        for i in range(len(c.cut_refs)):
            self.cut_selector.addItem(str(i + 1))
        self.cut_selector.currentIndexChanged.connect(self.select_cut)
        self.cut_selector.move(QPoint(265, 5))
        self.reverse_direction_checkbox = QCheckBox(parent=self, text='rev?')
        self.reverse_direction_checkbox.move(QPoint(300, 5))
        self.reverse_direction_checkbox.clicked.connect(self.reverse_direction)

        self.add_cut_button = QPushButton(parent=self, text='Add Cut')
        self.add_cut_button.move(QPoint(340, 5))
        self.add_cut_button.clicked.connect(self.add_cut)

        self.add_section_button = QPushButton(parent=self, text='+')
        self.add_section_button.setFixedSize(20, 20)
        self.add_section_button.move(QPoint(67, 0))
        self.add_section_button.clicked.connect(self.add_section)
        self.remove_section_button = QPushButton(parent=self, text='-')
        self.remove_section_button.setFixedSize(20, 20)
        self.remove_section_button.move(QPoint(88, 0))
        self.remove_section_button.clicked.connect(self.remove_section)

        self.set_section_positions()
        self.cut_selector.setCurrentIndex(0)

        self.show()

        self.reverse_direction_checkbox.setChecked(self.core.reverse_loaded[0])
        self.reverse_direction(self.reverse_direction_checkbox.isChecked())

    def set_section_positions(self):

        self.core.set_rect_positions(self.cut_selector.currentIndex())

        # t_h = 0.
        for i, r in enumerate(self.core.cs_rects):
            if i == 0:
                # r.moveTo(QPoint(50, 640 - r.height()))
                r.moveTo(QPoint(50, r.y()))
            else:
                dw = r.width() - self.core.cs_rects[0].width()
                # r.moveTo(QPoint(50 - int(dw / 2), 640 - t_h - r.height()))
                r.moveTo(QPoint(50 - int(dw / 2), r.y()))
            # t_h += r.height()

        self.force_refresh()

    def add_section(self):

        pml = self.core.project_manager_list

        if pml is None:
            last_sec = self.core.core_sections[-1]
            d = last_sec.diameter
            l = last_sec.length
            m = last_sec.mass
            self.core.add_section(d, l, m)
            self.core.set_rect_positions(self.cut_selector.currentIndex())
        else:
            pm = pml.parent()
            plw = pm.plug_list_widget
            l = plw.lists[0].count()
            if l:
                ps = PlugSelector(self)

                for i in range(l):
                    item = plw.lists[0].item(i)
                    item = QListWidgetItem(item.text())
                    ps.add_item(item)

                ps.show()

        self.force_refresh()

    def remove_section(self):

        self.core.remove_section()
        self.force_refresh()

    def reverse_direction(self, checked: bool):

        self.core.reverse_loaded[self.cut_selector.currentIndex()] = checked
        self.core.set_rect_positions(self.cut_selector.currentIndex())
        self.repaint()
        self.force_refresh()

    def add_cut(self):

        self.cut_selector.addItem(str(self.cut_selector.count() + 1))

        lengths = []
        for cs in self.core.core_sections:
            lengths.append(cs.length)
        self.core.create_cut(lengths, self.reverse_direction_checkbox.isChecked())

        self.cut_selector.setCurrentIndex(self.cut_selector.count() - 1)

    def select_cut(self):

        cut_num = self.cut_selector.currentIndex()
        for i, cs in enumerate(self.core.core_sections):
            cs.length = cs.cuts[self.core.cut_refs[cut_num][i]]
            r = self.core.cs_rects[i]
            self.set_core_entry_form(i, r)
            self.core_entry_form.send_edits()
        self.reverse_direction_checkbox.setChecked(self.core.reverse_loaded[cut_num])

    def force_refresh(self):

        w = self.width()
        h = self.height()
        self.setFixedSize(w + 1, h)
        self.setFixedSize(w, h)

    def collapse(self):
        self.setFixedSize(175, 670)
        self.form.size_to_widget()
        self.core_entry_form.current_section = -1
        self.core_entry_form.set_visible(False)
        self.add_section_button.setVisible(True)
        self.remove_section_button.setVisible(True)

    def paintEvent(self, a0: QPaintEvent):

        h_p = QPainter(self)
        h_p.setPen(QPen(Qt.black))
        h_p.setBrush(QBrush(QColor(200, 200, 200)))
        h_p.drawRect(self.core_holder.rect)

        p = QPainter(self)
        p.setBrush(QBrush(QColor(255, 248, 230), style=Qt.Dense3Pattern))
        p.setPen(QPen(Qt.black))

        for r in self.core.cs_rects:
            p.drawRect(r)

    def mousePressEvent(self, a0: QMouseEvent):

        x = a0.x()
        y = a0.y()

        for i, r in enumerate(self.core.cs_rects):
            if r.x() < x < (r.x() + r.width()) and r.y() < y < (r.y() + r.height()):

                self.set_core_entry_form(i, r)

    def set_core_entry_form(self, i: int, r: QRect):

        self.core_entry_form.current_section = i
        xe = self.core_holder.rect.width() + 40
        ye = r.y()
        self.setFixedWidth(xe + 170 + 155)
        self.core_entry_form.move(xe, ye)
        cs = self.core.core_sections[i]
        self.core_entry_form.length_edit.setText(str(cs.length))
        self.core_entry_form.diameter_edit.setText(str(cs.diameter))
        self.core_entry_form.mass_edit.setText(str(cs.mass))
        self.core_entry_form.porosity_edit.setText(str(cs.est_porosity))
        rho_string = '{:+.2f}'.format(cs.my_estimated_fluid_density())
        self.core_entry_form.fluid_density_display.setText(rho_string)
        self.core_entry_form.notes_edit.setText(cs.notes)

        option = 2
        if cs.lithology == 'sandstone':
            option = 0
        elif cs.lithology == 'limestone':
            option = 1

        self.core_entry_form.lithology_options.setCurrentIndex(option)

        self.core_entry_form.set_visible(True)
        self.add_section_button.setVisible(False)
        self.remove_section_button.setVisible(False)


class CoreEntryForm:

    def __init__(self, parent: CoreView):

        self.parent = parent

        self.length_label = QLabel('Length [cm]:', parent=parent)
        self.length_label.setFixedSize(100, 25)
        self.length_edit = QLineEdit(parent=parent)
        self.length_edit.setFixedSize(50, 25)
        self.length_edit.editingFinished.connect(self.send_edits)
        self.diameter_label = QLabel('Diameter [cm]:', parent=parent)
        self.diameter_label.setFixedSize(100, 25)
        self.diameter_edit = QLineEdit(parent=parent)
        self.diameter_edit.setFixedSize(50, 25)
        self.diameter_edit.editingFinished.connect(self.send_edits)
        self.diameter_edit.setEnabled(False)
        self.mass_label = QLabel('Init. Mass [g]:', parent=parent)
        self.mass_label.setFixedSize(100, 25)
        self.mass_edit = QLineEdit(parent=parent)
        self.mass_edit.setFixedSize(50, 25)
        self.mass_edit.editingFinished.connect(self.send_edits)
        self.mass_edit.setEnabled(False)
        self.porosity_label = QLabel(chr(981) + u' (est.)', parent=parent)
        self.porosity_label.setFixedSize(100, 25)
        self.porosity_edit = QLineEdit(parent=parent)
        self.porosity_edit.setFixedSize(50, 25)
        self.porosity_edit.editingFinished.connect(self.send_edits)
        self.porosity_edit.setEnabled(False)
        self.fluid_density_label = QLabel(chr(961) + u'_f (est.) [g/cc]:', parent=parent)
        self.fluid_density_label.setFixedSize(100, 25)
        self.fluid_density_display = QLabel('', parent=parent)
        self.fluid_density_display.setFixedSize(50, 25)
        self.lithology_label = QLabel(parent=parent, text='Lithology:')
        self.lithology_label.setFixedSize(100, 25)
        self.lithology_options = QComboBox(parent=parent)
        self.lithology_options.setFixedSize(75, 25)
        self.lithology_options.addItems(['Quartz', 'Calcite', 'Dolomite'])
        self.lithology_options.currentIndexChanged.connect(self.send_edits)
        self.lithology_options.setEnabled(False)
        self.notes_label = QLabel('Notes:', parent=parent)
        self.notes_label.setFixedSize(100, 25)
        self.notes_edit = JQTextEdit(parent=parent)
        self.notes_edit.setFixedSize(150, 175)
        self.notes_edit.focus_out_event.connect(self.send_edits)
        self.notes_edit.setEnabled(False)

        self.control_panel = [[self.length_label, self.diameter_label, self.mass_label,
                               self.porosity_label, self.fluid_density_label, self.lithology_label],
                              [self.length_edit, self.diameter_edit, self.mass_edit,
                               self.porosity_edit, self.fluid_density_display, self.lithology_options]]

        self.control_y = [y for y in range(0, 180, 30)]
        self.control_x = [[0, 0, 0, 0, 0, 0], [105, 105, 105, 105, 105, 80]]

        self.collapse_button = QPushButton(parent=parent, text='<<')
        self.collapse_button.setFixedSize(25, 25)
        self.collapse_button.clicked.connect(parent.collapse)

        # self.move(0, 0)

        self.current_section = -1

        self.set_visible(False)

    def set_visible(self, v: bool):

        self.length_label.setVisible(v)
        self.length_edit.setVisible(v)
        self.diameter_label.setVisible(v)
        self.diameter_edit.setVisible(v)
        self.mass_label.setVisible(v)
        self.mass_edit.setVisible(v)
        self.porosity_label.setVisible(v)
        self.porosity_edit.setVisible(v)
        self.fluid_density_label.setVisible(v)
        self.fluid_density_display.setVisible(v)
        self.lithology_label.setVisible(v)
        self.lithology_options.setVisible(v)
        self.notes_label.setVisible(v)
        self.notes_edit.setVisible(v)
        self.collapse_button.setVisible(v)

    def move(self, x: int, y: int):

        px = self.control_x
        py = self.control_y

        i = 0
        for ctl1, ctl2 in zip(self.control_panel[0], self.control_panel[1]):
            ctl1.move(QPoint(px[0][i] + x, py[i] + y))
            ctl2.move(QPoint(px[1][i] + x, py[i] + y))
            i += 1

        self.notes_label.move(QPoint(px[1][0] + self.control_panel[1][0].width() + 5 + x, y))
        self.notes_edit.move(QPoint(px[1][0] + self.control_panel[1][0].width() + 5 + x,
                                    self.notes_label.height() + 5 + y))

        self.collapse_button.move(px[1][4] + self.control_panel[1][0].width() -
                                  self.collapse_button.width() + x, py[-1] + py[1] - py[0] + y)

        if self.parent.height() < 200 + y:
            self.parent.setFixedHeight(225 + y)
            self.parent.form.size_to_widget()
        elif 200 + y <= 670:
            self.parent.setFixedHeight(670)
            self.parent.form.size_to_widget()

    def send_edits(self):

        if self.current_section == -1:
            return

        l = float(self.length_edit.text())
        d = float(self.diameter_edit.text())
        m = float(self.mass_edit.text())
        phi = float(self.porosity_edit.text())
        notes = self.notes_edit.toPlainText()
        i = self.lithology_options.currentIndex()
        lithology = 'sandstone'
        if i == 1:
            lithology = 'limestone'
        elif i == 2:
            lithology = 'dolomite'

        # cut_num = self.parent.cut_selector.currentIndex()

        self.parent.core.set_section_attributes(self.current_section, d, l, m, lithology, phi, notes)
        rho_string = '{:+.2f}'.format(self.parent.core.core_sections[self.current_section].my_estimated_fluid_density())
        self.fluid_density_display.setText(rho_string)
        self.parent.set_section_positions()

        r = self.parent.core.cs_rects[self.current_section]
        xe = self.parent.core_holder.rect.width() + 40
        ye = r.y()
        self.move(xe, ye)

        self.parent.force_refresh()


class JQTextEdit(QTextEdit):

    focus_out_event = pyqtSignal()

    def __init__(self, **kwargs):

        super(JQTextEdit, self).__init__(**kwargs)

    def focusOutEvent(self, e: QFocusEvent):

        super(JQTextEdit, self).focusOutEvent(e)
        self.focus_out_event.emit()


if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    form = CoreViewForm()
    sys.exit(app.exec_())
