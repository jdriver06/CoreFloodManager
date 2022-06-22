
from tkinter import Tk, Label

if __name__ == '__main__':
    app = Tk()
    my_label = Label(app, text='Hello world!')
    my_label.pack()
    app.mainloop()
