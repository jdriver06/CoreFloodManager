#!/usr/bin/env python

import sys
from os import path
# from datetime import datetime, timedelta
from win32com.client import Dispatch
from tkinter import Tk
import tkinter.messagebox as mbox


def print_to_dymo(strings: list, label_file: str):
    curdir = None
    if getattr(sys, 'frozen', False):
        # frozen
        curdir = path.dirname(sys.executable)
    else:
        # unfrozen
        curdir = path.dirname(path.abspath(__file__))

    mylabel = path.join(curdir, label_file)
    window = Tk()
    window.wm_withdraw()

    if not path.isfile(mylabel):
        mbox.showinfo('PyDymoLabel', 'Template file my.label does not exist')
        sys.exit(1)

    labelCom = Dispatch('Dymo.DymoAddIn')
    try:
        # now = datetime.now()
        # next = now + timedelta(30)
        labelCom = Dispatch('Dymo.DymoAddIn')
        labelText = Dispatch('Dymo.DymoLabels')
        isOpen = labelCom.Open(mylabel)
        selectPrinter = '\\\\LS300-PC\\DYMO LabelWriter 450'
        print(selectPrinter)
        labelCom.SelectPrinter(selectPrinter)

        # labelText.SetField('TEXT1', now.strftime('%Y/%m/%d'))
        # labelText.SetField('TEXT2', next.strftime('%Y/%m/%d'))
        i = 0
        for s in strings:
            # print('TEXT' + str(i + 1), strings[i])
            labelText.SetField('TEXT' + str(i + 1), strings[i])
            i += 1

        labelCom.StartPrintJob()
        labelCom.Print(1, False)
        labelCom.EndPrintJob()
    except:
        mbox.showinfo('PyDymoLabel', 'An error occurred during printing.')

    mbox.showinfo('PyDymoLabel', 'Label printed!')


if __name__ == "__main__":
    print_to_dymo(['hey', 'yo', 'this', 'test'], 'my.label')
