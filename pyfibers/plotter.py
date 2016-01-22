import os
class plotter(object):
    FIGURE_PATH = '/Users/troelsim/speciale/figures/'

    def __init__(self, name):
        self.name = name
        print "Plotting %s" % name

    def __enter__(self):
        self.plt = self.get_pyplotter()
        return self.plt

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.plt.savefig(os.path.join(self.FIGURE_PATH, self.name + self.extension()))
        pass

    def extension(self):
        return '.png'

    def get_pyplotter(self):
        pass


class pgfplotter(plotter):
    def get_pyplotter(self):
        import matplotlib
        matplotlib.use('pgf')
        matplotlib.rc('font', family='serif', size=12)
        matplotlib.rc('figure', figsize=[383.0/72, 327/72])
        from matplotlib import pyplot
        pyplot.clf()
        return pyplot

    def extension(self):
        return '.pgf'



if __name__=='__main__':
    with pgfplotter('testplot') as plt:
        plt.plot([1,2,3],[3,4,5])
