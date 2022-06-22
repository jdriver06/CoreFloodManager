
from sys import argv, exit
from os import listdir
from os.path import isfile, exists
from PyQt5.Qt import QApplication
from PyQt5.QtWidgets import QDialog, QFileDialog, QPushButton, QGridLayout, QLabel


__version__ = '1.0.0'


class InstallTool(QDialog):

    def __init__(self):
        super(InstallTool, self).__init__()
        self.setWindowTitle('Install Tool v' + __version__)

        root_dir = ''
        if exists('U:\\UEORLAB1\\Python Files\\ProjectManager\\'):
            root_dir = 'U:\\UEORLAB1\\Python Files\\ProjectManager\\'
        elif exists('Q:\\UEORLAB1\\Python Files\\ProjectManager\\'):
            root_dir = 'Q:\\UEORLAB1\\Python Files\\ProjectManager\\'

        self.root_dir = root_dir
        self.source_dir = ''
        self.destination_dir = ''

        self.setFixedSize(500, 150)
        lyt = QGridLayout()

        self.source_dir_btn = QPushButton(parent=self, text='Source')
        self.source_dir_btn.clicked.connect(self.source_btn_clicked)
        self.source_text = QLabel(parent=self, text='')
        self.destination_dir_btn = QPushButton(parent=self, text='Destination')
        self.destination_dir_btn.clicked.connect(self.destination_btn_clicked)
        self.destination_text = QLabel(parent=self, text='')
        self.install_btn = QPushButton(parent=self, text='Install')
        self.install_btn.clicked.connect(self.install_btn_clicked)

        lyt.addWidget(self.source_dir_btn, 0, 0, 1, 1)
        lyt.addWidget(self.destination_dir_btn, 1, 0, 1, 1)
        lyt.addWidget(self.install_btn, 2, 0, 1, 1)
        lyt.addWidget(self.source_text, 0, 1, 1, 2)
        lyt.addWidget(self.destination_text, 1, 1, 1, 2)

        self.setLayout(lyt)
        self.show()

    def update_displayed_text(self):

        self.source_text.setText(self.source_dir)
        self.destination_text.setText(self.destination_dir)

    def source_btn_clicked(self):

        self.source_dir = self.get_directory()
        self.update_displayed_text()

    def destination_btn_clicked(self):

        self.destination_dir = self.get_directory()
        self.update_displayed_text()

    def get_directory(self) -> str:

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        return QFileDialog.getExistingDirectory(self, directory=self.root_dir)

    def install_btn_clicked(self):

        directory_dict = InstallTool.get_directory_dict(self.source_dir)
        print(directory_dict)

    @staticmethod
    def get_directory_dict(directory: str) -> dict:

        directory_dict = {}

        result = listdir(directory)
        items = []

        for item in result:

            if isfile(directory + '\\' + item):
                items.append(item)
            else:
                items.append(InstallTool.get_directory_dict(directory + '\\' + item))

        directory_dict[directory] = items

        return directory_dict


def main():
    app = QApplication(argv)
    tool = InstallTool()
    print(tool)
    exit(app.exec_())


if __name__ == '__main__':
    main()
