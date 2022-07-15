
from produced_fluids import AqueousTracer
import numpy as np


if __name__ == '__main__':

    volumes = np.array([4.050, 5.100, 5.120, 5.500, 5.650, 6.000, 7.100, 7.050, 6.500, 6.700, 6.700, 7.050, 6.95, 6.90,
                        6.90, 6.80, 7.80, 6.80, 7.10, 7.30, 6.80, 6.80, 7.05, 6.90, 7.00, 7.00, 6.85, 6.70, 6.75, 6.92,
                        6.75, 7.00, 6.85, 7.00, 6.75, 6.75, 6.85, 6.72, 6.72, 7.40, 6.72, 6.60, 6.80, 7.00, 6.77, 7.00,
                        7.10, 6.90])

    aq_volumes = np.array([2.75, 5.10, 5.10, 5.50, 5.65, 6.00, 6.80, 6.80, 6.40, 6.700, 6.700, 6.95, 6.85, 6.85, 6.60,
                           6.05, 6.60, 6.20, 6.85, 7.10, 6.00, 6.60, 6.70, 6.80, 6.50, 6.65, 6.65, 6.70, 6.70, 6.90,
                           6.70, 6.40, 6.80, 6.70, 6.40, 6.70, 6.80, 6.70, 6.70, 7.00, 6.70, 6.60, 6.80, 7.00, 6.75,
                           6.75, 6.90, 6.65])

    pHs = np.array([np.nan, np.nan, 8.45, np.nan, 8.39, np.nan, 8.25, np.nan, 7.99, np.nan, 8.07, np.nan, 8.23, np.nan,
                    9.14, np.nan, 10.01, np.nan, 10.19, np.nan, 10.14, np.nan, 10.14, np.nan, 10.16, np.nan, 9.91,
                    np.nan, 9.89, np.nan, 9.90, np.nan, 9.68, np.nan, 9.78, np.nan, 9.61, np.nan, 9.57, np.nan, 9.65,
                    np.nan, 9.60, np.nan, np.nan, np.nan, np.nan, np.nan])

    OH_conc = np.power(10., pHs - 14.)
    no_nans = np.where(~np.isnan(OH_conc))

    tracer = AqueousTracer('Alkali Consumption', np.cumsum(volumes), np.cumsum(aq_volumes), OH_conc)
    tracer.set_normalization_max(1.)

    recovered = 0.001 * tracer.get_area_above_below_tracer(False)

    print('[OH-]: {:.2E} excess moles recovered.'.format(recovered))

    injected = 0.001 * 46.1 * np.power(10., 10.62 - 14.)

    print('injected: {:.2E}'.format(injected))

    print('consumed: {:.1f}%'.format(100. * (injected - recovered) / injected))
