
from sys import argv, exit
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QPushButton, QLabel, QFileDialog, QSpinBox, QHBoxLayout
from PyQt5.QtChart import QChart, QChartView, QValueAxis, QLineSeries
from PyQt5.Qt import QApplication, Qt, QColor, QIcon, QPointF
import csv


__version__ = '1.1.0'


class PressureDataAveragerGUI(QMainWindow):

    def __init__(self):
        super(PressureDataAveragerGUI, self).__init__()
        self.setFixedSize(600, 600)

        self.setWindowTitle('Pressure Data Averager')
        self.setWindowIcon(QIcon('UEORS logo cropped.png'))
        self.file_path = ''
        self.period = -1

        self.c_widget = PressureDataAveragerWidget(self)
        self.setCentralWidget(self.c_widget)

        self.p_data = [[], [], [], [], []]
        self.average_data = [[], [], [], [], []]
        self.t_data = []

        self.show()


class PressureDataAveragerWidget(QWidget):

    def __init__(self, parent: PressureDataAveragerGUI):
        super(PressureDataAveragerWidget, self).__init__(parent=parent)
        self.pd_gui = parent

        self.setLayout(QVBoxLayout())
        self.layout().setSpacing(0)

        self.open_button = QPushButton(parent=self, text='Open')
        font = self.open_button.font()
        font.setPointSize(11)
        self.open_button.setFont(font)
        self.open_button.clicked.connect(self.open_file)

        self.save_button = QPushButton(parent=self, text='Save As')
        self.save_button.setFont(font)
        self.save_button.clicked.connect(self.write_file)

        self.file_text = QLabel(parent=self, text='')
        self.file_text.setFont(font)

        self.period_text = QLabel(parent=self, text='')
        self.period_text.setFont(font)
        self.points_spin_label = QLabel(parent=self, text='Points to Average: ')
        self.points_spin_label.setFont(font)
        self.points_spin = QSpinBox(parent=self)
        self.points_spin.setMinimum(1)
        self.points_spin.setFont(font)
        self.recalculate_button = QPushButton(parent=self, text=' Recalculate ')
        self.recalculate_button.setFont(font)
        self.recalculate_button.clicked.connect(self.display_data)
        self.chart = QChart(flags=Qt.WindowFlags())
        self.chart.legend().hide()
        self.x_axis = QValueAxis()
        self.y_axis = QValueAxis()

        b_widget = QWidget(parent=self)
        b_widget.setLayout(QHBoxLayout())
        b_widget.layout().addWidget(self.open_button)
        b_widget.layout().addWidget(self.save_button)

        f_widget = QWidget(parent=self)
        f_widget.setLayout(QHBoxLayout())
        f_widget.layout().addWidget(self.file_text)

        p_widget = QWidget(parent=self)
        p_widget.setLayout(QHBoxLayout())
        p_widget.layout().addWidget(self.period_text)

        a_widget = QWidget(parent=self)
        a_widget.setLayout(QHBoxLayout())
        a_widget.layout().addWidget(self.points_spin_label)
        a_widget.layout().addWidget(self.points_spin)
        a_widget.layout().addWidget(self.recalculate_button)
        a_widget.setFixedWidth(280)

        self.curves = [QLineSeries(), QLineSeries(), QLineSeries(), QLineSeries(), QLineSeries()]
        for curve in self.curves:
            self.chart.addSeries(curve)
            curve.attachAxis(self.x_axis)
            curve.attachAxis(self.y_axis)

        self.chart.setAxisX(self.x_axis)
        self.chart.setAxisY(self.y_axis)
        self.view = QChartView(self.chart)

        self.layout().addWidget(b_widget)
        self.layout().addWidget(f_widget)
        self.layout().addWidget(p_widget)
        self.layout().addWidget(a_widget)
        self.layout().addWidget(self.view)

        self.update_file_info()

    def open_file(self):

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_path, _ = QFileDialog.getOpenFileName(self.parent(),
                                                   'Open file...',
                                                   'U:\\UEORLAB1\\',
                                                   'TXT Files (*.txt)',
                                                   options=options)

        if not file_path:
            return

        self.pd_gui.file_path = file_path
        file_path_split = file_path.split('/')
        self.file_text.setText(file_path_split[-1])

        w_data = []
        s1_data = []
        s2_data = []
        s3_data = []
        s4_data = []

        period = 1.

        with open(file_path, 'r') as f:

            reader = csv.reader(f, delimiter='\t')
            eof = False
            i = 0
            while not eof:
                try:
                    result = next(reader)
                    if i == 2:
                        period_text = result[0]
                        period = float(str.split(period_text)[-1])
                    if i > 11:
                        w_data.append(float(result[0]))
                        s1_data.append(float(result[1]))
                        s2_data.append(float(result[2]))
                        s3_data.append(float(result[3]))
                        s4_data.append(float(result[4]))
                except Exception as e:
                    print(e)
                    if i > 11:
                        eof = True
                i += 1

        self.pd_gui.period = period
        t_data = [period * i / 60.0 for i in range(len(w_data))]

        self.pd_gui.t_data = t_data
        self.pd_gui.p_data = [w_data, s1_data, s2_data, s3_data, s4_data]

        self.points_spin.setValue(1)
        self.display_data()

    def display_data(self):

        if not self.pd_gui.p_data[0]:
            return

        if len(self.pd_gui.p_data[0]) < self.points_spin.value():
            return

        colors = ['#4F81BD', '#C0504D', '#9BBB59', '#8064A2', '#F79646']

        y_min = self.y_axis.min()
        y_max = self.y_axis.max()

        n = self.points_spin.value()
        average_data = []
        offset = self.pd_gui.period * n / 120.0

        for curve, data, color in zip(self.curves, self.pd_gui.p_data, colors):

            a_data = []
            for i in range(0, len(data), n):
                if i + n <= len(data):
                    a_data.append(sum(data[i:i + n]) / n)

            average_data.append(a_data)
            t_data = [self.pd_gui.period * n * j / 60.0 for j in range(len(a_data))]

            d_min = min(a_data)
            d_max = max(a_data)
            if d_min < y_min:
                y_min = d_min
            if d_max > y_max:
                y_max = d_max

            curve.clear()
            k = 0
            for t, a in zip(t_data, a_data):
                curve.append(QPointF(t + offset, a))
                k += 1

            curve.attachAxis(self.x_axis)
            curve.attachAxis(self.y_axis)

            pen = curve.pen()
            pen.setColor(QColor(color))

            curve.setPen(pen)

        self.pd_gui.average_data = average_data

        self.x_axis.setMin(0.)
        self.x_axis.setMax(max(self.pd_gui.t_data) + offset)
        self.y_axis.setMin(y_min)
        self.y_axis.setMax(y_max)

        self.x_axis.setTitleText('t [min]')
        self.y_axis.setTitleText(chr(916) + u'P [psi]')

        self.update_file_info()

    def update_file_info(self):

        file_path_split = self.pd_gui.file_path.split('/')
        self.file_text.setText('File: {}'.format(file_path_split[-1]))

        if self.pd_gui.period == -1:
            period = ''
            new_period = ''
        else:
            period = str(self.pd_gui.period)
            new_period = str(self.pd_gui.period * self.points_spin.value())

        self.period_text.setText('Period: {}s, New Period: {}s'.format(period, new_period))

    def write_file(self):

        if not self.pd_gui.file_path:
            return

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_path, _ = QFileDialog.getSaveFileName(self.parent(), 'Save as...', '', 'TXT File (*.txt)', options=options)

        if not file_path:
            return

        if not self.pd_gui.p_data[0]:
            return

        data = self.pd_gui.average_data

        w = data[0]
        s1 = data[1]
        s2 = data[2]
        s3 = data[3]
        s4 = data[4]

        with open(file_path, 'w') as f:

            writer = csv.writer(f, delimiter='\t')

            for c0, c1, c2, c3, c4 in zip(w, s1, s2, s3, s4):
                writer.writerow([c0, c1, c2, c3, c4])


if __name__ == '__main__':
    app = QApplication(argv)
    gui = PressureDataAveragerGUI()
    exit(app.exec_())
