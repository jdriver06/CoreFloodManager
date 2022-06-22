
from reportlab.pdfgen import canvas
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter, portrait, landscape
from reportlab.lib.units import inch
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Table, TableStyle, KeepTogether


def plan_first_page(can: canvas.Canvas, doc: SimpleDocTemplate):

    can.saveState()

    can.setFont('Helvetica-Bold', 14)
    title_str = doc.filename[:-4]
    title_str = title_str.split('/')
    title_str = title_str[-1]
    can.drawCentredString(4.25 * inch, 10.1 * inch, title_str)

    doc.afterPage()
    can.restoreState()


def plan_later_pages(can: canvas.Canvas, _):

    can.saveState()
    can.setFont('Helvetica', 12)
    can.restoreState()


def create_plan(fe, file_name: str):

    doc = SimpleDocTemplate(file_name + '.pdf',
                            leftMargin=0.75*inch,
                            righMargin=0.75*inch,
                            topMargin=1.*inch,
                            bottomMargin=1.*inch)

    doc.pagesize = portrait(letter)

    story = [Paragraph('Core to Use: ' + fe.core.name + '<br />\n'),
             Paragraph('Objective: ' + fe.objective + '<br />\n<br />\n')]

    for fl in fe.floods:

        flood_text = ''

        flood_text += '<b>' + fl.name + '</b><br />\n'
        flood_text += 'Temperature [C]: {:.1f}, Back Pressure [psig]: {:.1f}, Planned Volume [mL]: {:.1f}'.\
            format(fl.temperature, fl.back_pressure, fl.planned_volume)
        flood_text += '<br />\n'
        if fl.notes:
            flood_text += 'Plan Notes: ' + fl.notes + '<br />\n'
        if fl.fluid is not None:
            flood_text += 'Injection Fluid: ' + fl.fluid.name + '<br />\n'
        else:
            flood_text += 'Injection Fluid: None Selected<br />\n'
        flood_text += '<br />\n'

        story.append(KeepTogether(Paragraph(flood_text)))

    fluid_text = '<br />\n<b>Injection Fluids:</b><br />\n'

    for flu in fe.injection_fluids_list.signal_lists[0].objects:

        fluid_text += flu.name + ' ('

        if hasattr(flu.specific_fluid, 'brine'):
            fluid_text += flu.specific_fluid.brine.name
        else:
            oil_sample = flu.specific_fluid.oil_sample
            fluid_text += oil_sample.ref_objects[0].name
            for additive, c in zip(oil_sample.additives, oil_sample.concentrations):
                fluid_text += ', {}: {} ppm'.format(additive.name, int(c))

        for additive, c in zip(flu.specific_fluid.additives, flu.concentrations):
            fluid_text += ', {}: {} ppm'.format(additive.name, int(c))

        fluid_text += ')<br />\n'

    story.append(KeepTogether(Paragraph(fluid_text)))

    doc.build(story, onFirstPage=plan_first_page, onLaterPages=plan_later_pages)


def create_injection_fluid_mixing_sheet(stock_names_list: list, stocks_mixing_list: list, balance_brine: list,
                                        fluid_name: str, mixing_list: list, salinity: float, file_name: str):

    doc = SimpleDocTemplate(file_name + '.pdf',
                            leftMargin=1.5 * inch,
                            righMargin=1.5 * inch,
                            topMargin=1. * inch,
                            bottomMargin=1. * inch)

    doc.pagesize = landscape(letter)

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name='centered', alignment=TA_CENTER))

    story = []
    stock_tables = []
    stock_pH_tables = []
    stock_titles = []

    if stocks_mixing_list:
        for stock_name, stock_mixing in zip(stock_names_list, stocks_mixing_list):
            stock_titles.append(Paragraph('<b>' + stock_name + '</b>', style=styles['centered']))
            stocks_data = [['Component', 'Lot Number', 'Stock [%]', 'Target [ppm]', 'Target [g]', 'Actual [g]']]
            total_mass = 0.
            for item in stock_mixing:
                component = item[0]
                lot_num = item[1]
                stock_c = item[2]
                if stock_c == 1000000.:
                    stock_c = '100.0'
                else:
                    stock_c = '{:.1f}'.format(1.e-4 * stock_c)
                target_c = item[3]
                if target_c == 1000000.:
                    target_c = ''
                else:
                    target_c = '{:.0f}'.format(target_c)
                target_mass = '{:.3f}'.format(item[4])
                total_mass += item[4]
                actual_mass = ''
                stocks_data.append([component, lot_num, stock_c, target_c, target_mass, actual_mass])

            stocks_data.append(['Total', '', '', '', '{:.3f}'.format(total_mass), ''])
            st = Table(data=stocks_data)
            st.setStyle(TableStyle([('INNERGRID', (0, 1), (-1, -2), 0.25, colors.black),
                                    ('INNERGRID', (4, -1), (4, -1), 0.25, colors.black),
                                    ('BOX', (0, 1), (-1, -2), 0.25, colors.black),
                                    ('BOX', (4, -1), (4, -1), 0.25, colors.black),
                                    ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
                                    ('FONTNAME', (0, 0), (-1, 0), 'Courier-Bold'),
                                    ('FONTNAME', (0, -1), (-1, -1), 'Courier-Bold'),
                                    ('FONTNAME', (0, 1), (-1, -2), 'Courier')]))

            stock_tables.append(st)
            data_pH = [['', ''], ['pH Init:', '__________'], ['pH Final:', '__________'], ['', '']]
            pH = Table(data=data_pH, hAlign='RIGHT')
            pH.setStyle(TableStyle([('ALIGN', (0, 0), (-1, -1), 'RIGHT')]))
            stock_pH_tables.append(pH)

    bbt = None
    if balance_brine:
        bb_data = []
        ordered_keys = ["Li+", "Na+", "K+", "Mg++", "Ca++", "Ba++", "Sr++", "Fe++", "F-", "Cl-", "Br-", "HCO3-",
                        "CO3--", "SO4--"]
        for key, val in zip(ordered_keys, balance_brine):
            bb_data.append([key, round(val, 1)])

        bbt = Table(data=bb_data)
        bbt.setStyle(TableStyle([('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black),
                                 ('BOX', (0, 0), (-1, -1), 0.25, colors.black),
                                 ('FONTNAME', (0, 0), (0, -1), 'Courier-Bold'),
                                 ('FONTNAME', (1, 0), (1, -1), 'Courier')]))

    data = [['Component', 'Lot Number', 'Stock [ppm]', 'Target [ppm]', 'Target [g]', 'Actual [g]']]
    total_mass = 0.
    for item in mixing_list:
        component = item[0]
        lot_num = item[1]
        stock_c = item[2]
        if stock_c == 1000000.:
            stock_c = 'pure'
        else:
            stock_c = '{:.0f}'.format(stock_c)
        target_c = item[3]
        if target_c == 1000000.:
            target_c = ''
        else:
            target_c = '{:.0f}'.format(target_c)
        target_mass = '{:.3f}'.format(item[4])
        total_mass += item[4]
        actual_mass = ''
        data.append([component, lot_num, stock_c, target_c, target_mass, actual_mass])
    data.append(['Total (Salinity = {:.0f} TDS)'.format(salinity), '', '', '', '{:.3f}'.format(total_mass), ''])

    t = Table(data=data)
    t.setStyle(TableStyle([('INNERGRID', (0, 1), (-1, -2), 0.25, colors.black),
                           ('INNERGRID', (4, -1), (4, -1), 0.25, colors.black),
                           ('BOX', (0, 1), (-1, -2), 0.25, colors.black),
                           ('BOX', (4, -1), (4, -1), 0.25, colors.black),
                           ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
                           ('FONTNAME', (0, 0), (-1, 0), 'Courier-Bold'),
                           ('FONTNAME', (0, -1), (-1, -1), 'Courier-Bold'),
                           ('FONTNAME', (0, 1), (-1, -2), 'Courier')]))

    data_pHRI = [['', ''], ['pH Init:', '__________'], ['pH Final:', '__________'], ['R.I. [o/oo]:', '__________'],
                 ['Visc. [cP] @ ___ s-1', '__________']]
    pHRI = Table(data=data_pHRI, hAlign='RIGHT')
    pHRI.setStyle(TableStyle([('ALIGN', (0, 0), (-1, -1), 'RIGHT')]))

    for title, table, pH_table in zip(stock_titles, stock_tables, stock_pH_tables):
        stock_fluid = [title, table, pH_table]
        story.append(KeepTogether(stock_fluid))
    if bbt is not None:
        story.append(KeepTogether([Paragraph('Balance Brine', style=styles['centered']), bbt]))
    last_fluid = [Paragraph('<b>' + fluid_name + '</b>', style=styles['centered']), t, pHRI]
    story.append(KeepTogether(last_fluid))

    doc.build(story)
