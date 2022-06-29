
from sys import argv, exit
from PyQt5.Qt import QApplication, QIcon, QThread, QObject, pyqtSignal
from PyQt5.QtWidgets import QMainWindow, QWidget, QGridLayout, QLabel, QPushButton, QTextEdit
from os.path import isdir, isfile
from os import listdir, getcwd, startfile, remove
import shutil
from time import sleep
import db_tools

__version__ = '1.1.1'


def is_version_greater(v1: str, v2: str) -> bool:

    v1 = v1.split('.')
    v2 = v2.split('.')

    for i in range(3):
        if int(v2[i]) > int(v1[i]):
            return False
        elif int(v1[i]) > int(v2[i]):
            return True

    return False


def bubble_sort_directory_list_by_version(directories: list) -> list:

    is_sorted = False

    while not is_sorted:
        is_sorted = True

        for i in range(len(directories) - 1):
            di = directories[i]
            dj = directories[i + 1]
            if is_version_greater(di.split('v')[-1].replace('-', '.'), dj.split('v')[-1].replace('-', '.')):
                directories[i] = dj
                directories[i + 1] = di
                is_sorted = False

    return directories


class UpdateWorker(QObject):

    started = pyqtSignal()
    finished = pyqtSignal()
    progress = pyqtSignal(int, float)

    def __init__(self, u_or_q: str, upgrade_path: list, gui: str):
        super(UpdateWorker, self).__init__()
        self.upgrade_path = upgrade_path
        self.u_or_q = u_or_q
        self.gui = gui

    def update_to_current(self):

        sleep(0.1)
        self.started.emit()

        path = self.u_or_q + ':\\UEORLAB1\\Python Files\\' + self.gui + '\\'
        full_upgrade_path = [path + self.gui + '_v' + up.replace('.', '-') for up in self.upgrade_path]
        cwd = getcwd()
        current_files_and_directories = listdir(cwd)

        for i, fup in enumerate(full_upgrade_path):
            all_files_and_directories = listdir(fup)
            l = len(all_files_and_directories)
            for j, f_or_d in enumerate(all_files_and_directories):
                self.progress.emit(i + 1, j / l)
                if f_or_d not in current_files_and_directories:
                    sleep(0.01)
                    if isfile(fup + '\\' + f_or_d) and not isfile(cwd + '\\' + f_or_d):
                        shutil.copy2(fup + '\\' + f_or_d, cwd)
                    elif isdir(fup + '\\' + f_or_d) and not isdir(cwd + '\\' + f_or_d):
                        shutil.copytree(fup + '\\' + f_or_d, cwd + '\\' + f_or_d)

        sleep(0.1)
        self.finished.emit()


class GUI_Launcher(QMainWindow):

    def __init__(self, cnxn: db_tools.pyodbc.Connection, gui: str, icon_file: str):
        super(GUI_Launcher, self).__init__()
        self.cnxn = cnxn
        self.gui = gui
        self.icon_file = icon_file
        self.setFixedSize(600, 200)
        self.setWindowTitle('{} Launcher'.format(gui))
        self.setWindowIcon(QIcon(icon_file))

        self.upgrade_widget = QWidget()
        lyt = QGridLayout()
        self.upgrade_widget.setLayout(lyt)
        self.upgrade_label = QLabel(parent=self.upgrade_widget, text='Software up-to-date.')
        font = self.upgrade_label.font()
        font.setPointSize(11)
        self.upgrade_label.setFont(font)
        self.upgrade_button = QPushButton(parent=self.upgrade_widget, text='Update')
        self.upgrade_button.setFont(font)
        self.upgrade_button.clicked.connect(self.update_clicked)
        self.upgrade_button.setEnabled(False)
        self.upgrade_button.setFixedWidth(100)
        self.launch_button = QPushButton(parent=self.upgrade_widget, text='Launch')
        self.launch_button.setFont(font)
        self.launch_button.clicked.connect(self.launch)
        self.launch_button.setEnabled(False)
        self.launch_button.setFixedWidth(100)

        self.upgrade_text = QTextEdit(parent=self)
        self.upgrade_text.setEnabled(False)
        self.left_button = QPushButton(parent=self, text='<')
        self.left_button.setFont(font)
        self.left_button.clicked.connect(self.left_button_clicked)
        self.left_button.setEnabled(False)
        self.right_button = QPushButton(parent=self, text='>')
        self.right_button.setFont(font)
        self.right_button.clicked.connect(self.right_button_clicked)
        self.right_button.setEnabled(False)

        self.upgrade_widget.setContentsMargins(0, 0, 0, 0)
        lyt.addWidget(self.upgrade_label, 0, 0, 1, 1)
        lyt.addWidget(self.upgrade_button, 0, 1, 1, 1)
        lyt.addWidget(self.launch_button, 0, 2, 1, 1)
        lyt.addWidget(self.upgrade_text, 1, 0, 1, 3)
        lyt.addWidget(self.left_button, 2, 1, 1, 1)
        lyt.addWidget(self.right_button, 2, 2, 1, 1)

        self.setCentralWidget(self.upgrade_widget)

        self.u_or_q = 'U'
        if not isdir('U:\\'):
            self.u_or_q = 'Q'

        self.upgrade_path = []
        self.file = ''
        self.info = []
        self.i = 0
        self.worker = None
        self.t = None

    def check_most_recent_exe(self) -> tuple:

        if not self.file:
            return ()

        ver = self.file.split('v')[-1].split('.exe')[0].split('-')

        DBO = db_tools.DBObjects
        result = db_tools.execute_stored_procedure_with_params(self.cnxn, DBO.Procedures.CheckPMUpdateAvailable.value,
                                                               {DBO.Params.v1.value: ver[0],
                                                                DBO.Params.v2.value: ver[1],
                                                                DBO.Params.v3.value: ver[2]}).fetchall()[0]

        return result

    def find_most_recent_exe(self):

        my_dir = getcwd()
        results = listdir(my_dir)
        exe = ''

        for result in results:

            if result.find(self.gui + '_v') == -1 or result.find('.manifest') != -1:
                continue

            if isfile(my_dir + '\\' + result) and result > exe:
                exe = result

        self.file = my_dir + '\\' + exe

        new_exe = self.check_most_recent_exe()
        old_exe = exe
        exe = exe.split('v')[-1].split('.')[0].split('-')
        exe = (int(exe[0]), int(exe[1]), int(exe[2]))

        if new_exe[0] is not None:
            upgrade_path = self.check_for_update(old_exe)
            self.upgrade_path = upgrade_path
            self.upgrade_label.setText('Software needs updating to version {}.{}.{}.'.format(*new_exe))
            self.upgrade_label.setStyleSheet('QLabel {color: red};')
            self.upgrade_button.setEnabled(True)
            DBO = db_tools.DBObjects
            info = db_tools.execute_stored_procedure_with_params(self.cnxn, DBO.Procedures.GetPMVerNotesCurr.value,
                                                                 {DBO.Params.v1.value: new_exe[0],
                                                                  DBO.Params.v2.value: new_exe[1],
                                                                  DBO.Params.v3.value: new_exe[2],
                                                                  DBO.Params.cv1.value: exe[0],
                                                                  DBO.Params.cv2.value: exe[1],
                                                                  DBO.Params.cv3.value: exe[2]}).fetchall()

            self.info = info

            is_required = False
            for piece in info:
                if piece[4]:
                    is_required = True
                    break

            self.launch_button.setEnabled(not is_required)

            self.i = len(info) - 1
            self.load_info()
            self.right_button.setEnabled(True)
            self.left_button.setEnabled(True)

        else:
            self.launch_button.setEnabled(True)
            self.right_button.setEnabled(False)
            self.left_button.setEnabled(False)

    def right_button_clicked(self):

        if self.i == len(self.info) - 1:
            return

        self.i += 1
        self.load_info()

    def left_button_clicked(self):

        if self.i == 0:
            return

        self.i -= 1
        self.load_info()

    def load_info(self):

        piece = self.info[self.i]
        required = 'Yes'
        if not piece[4]:
            required = 'No'
        text = 'Version: {}.{}.{}\tRequired? {}.\n\nNotes: {}'
        self.upgrade_text.setText(text.format(piece[0], piece[1], piece[2], required, piece[3]))

    def remove_previous_exes(self):

        my_dir = getcwd()
        results = listdir(my_dir)
        exe = self.file.split('\\')[-1].split('.')[0]

        for result in results:

            if result.find(self.gui + '_v') > -1 and result.find(exe) == -1:
                remove(my_dir + '\\' + result)

    def check_for_update(self, exe) -> list:

        v = exe.split('v')[-1]
        v = v.split('.')[0]
        v = v.replace('-', '.')

        path = self.u_or_q + ':\\UEORLAB1\\Python Files\\' + self.gui + '\\'
        directories = [f for f in listdir(path) if not isfile(path + f) and f.find(self.gui + '_v') != -1]
        directories = bubble_sort_directory_list_by_version(directories)

        latest_version = '0.0.0'
        upgrade_path = []

        for d in directories:
            version = d.split('v')[-1]
            version = version.replace('-', '.')

            if is_version_greater(version, latest_version):
                latest_version = version
                if is_version_greater(version, v):
                    upgrade_path.append(version)

        return upgrade_path

    def update_clicked(self):

        self.launch_button.setEnabled(False)
        self.worker = UpdateWorker(self.u_or_q, self.upgrade_path, self.gui)
        self.t = QThread()
        self.worker.moveToThread(self.t)
        self.t.started.connect(self.worker.update_to_current)
        self.worker.finished.connect(self.t.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.worker.finished.connect(self.update_done)
        self.worker.progress.connect(self.update_progress)
        self.worker.started.connect(self.update_started)
        self.t.finished.connect(self.t.deleteLater)
        self.t.start(priority=QThread.LowPriority)

    def update_started(self):

        self.upgrade_button.setEnabled(False)
        self.upgrade_label.setText('Software updating. Do not close.')

    def update_progress(self, i: int, p: float):

        self.upgrade_label.setText('Update progress: {} of {}, {:.1f}%'.format(i, len(self.upgrade_path), 100. * p))

    def update_done(self):

        self.find_most_recent_exe()
        self.upgrade_label.setText('Removing previous exe files.')
        self.remove_previous_exes()
        self.upgrade_label.setText('Software ready to launch.')
        self.upgrade_label.setStyleSheet('QLabel {color: black};')
        self.launch_button.setEnabled(True)
        self.upgrade_button.setEnabled(False)

    def launch(self):

        self.find_most_recent_exe()
        startfile(self.file)
        self.close()

    def close(self):
        self.cnxn.close()
        super(GUI_Launcher, self).close()


def main():

    server = 'UEORS-DB\\'
    database = 'UEORS_MAIN_DB'
    username = 'ueors_user'
    password = 'ueor101'

    try:
        cnxn = db_tools.connect_to_database(server, database, username, password)
    except Exception as e:
        print('error with cnxn.')
        print(e)
        return

    try:
        app = QApplication(argv)
        app.setStyle('fusion')
        gui = 'ProjectManager'
        icon_file = 'Coreholder Cropped.jpg'
        launcher = GUI_Launcher(cnxn=cnxn, gui=gui, icon_file=icon_file)
        launcher.show()
        launcher.find_most_recent_exe()
        exit(app.exec_())
    except Exception as e:
        print('error with launcher init.')
        print(e)
        return


if __name__ == '__main__':
    main()
