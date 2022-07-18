
import j_utils


class Dummy(j_utils.SignalView):

    def __init__(self):
        super(Dummy, self).__init__()
        self.attr = 3


if __name__ == '__main__':
    print(hasattr(Dummy(), 'close_signal'))
