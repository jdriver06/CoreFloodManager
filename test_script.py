
from PyQt5.Qt import QApplication
from sys import argv, exit


if __name__ == '__main__':

    hey = {'good': 5}
    hey.pop('good')
    print(hey)

    app = QApplication(argv)
    psc = app.primaryScreen()
    dpi = psc.physicalDotsPerInch()
    dpi_x = psc.physicalDotsPerInchX()
    dpi_y = psc.physicalDotsPerInchY()
    p_size = psc.physicalSize()
    geo = psc.geometry()
    print('dpi: {}, dpi_x: {}, dpi_y: {}, physical size: {}, geometry: {}'.format(dpi, dpi_x, dpi_y, p_size, geo))
    exit(app.exit())
