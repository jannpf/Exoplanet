import sys


class Data:
    _instance = None

    def __init__(self):
        self.t = []
        self.y = []
        self.sig = []

    def load(self, filename):
        try:
            with open(filename, "r") as fin:
                lines = fin.readlines()
                for line in lines:
                    temp1, temp2, temp3 = map(float, line.split())
                    self.t.append(temp1)
                    self.y.append(temp2)
                    self.sig.append(temp3)
            print(f"# Loaded {len(self.t)} data points from file {filename}.")
        except IOError:
            print(f"# Error. Couldn't open file {filename}.")
            sys.exit(1)

    # Getters
    def get_t(self):
        return self.t

    def get_y(self):
        return self.y

    def get_sig(self):
        return self.sig

    def get_y_min(self):
        return min(self.y)

    def get_y_max(self):
        return max(self.y)

    # Singleton instance
    @staticmethod
    def get_instance():
        if Data._instance is None:
            Data._instance = Data()
        return Data._instance


if __name__ == "__main__":
    data = Data.get_instance()
    data.load("data/nu_oph.txt")  # Replace with your actual data file name

    t_data = data.get_t()
    y_data = data.get_y()
    sig_data = data.get_sig()
    y_min = data.get_y_min()
    y_max = data.get_y_max()

    print(f"t data: {t_data}")
    print(f"y data: {y_data}")
    print(f"sig data: {sig_data}")
    print(f"Minimum y: {y_min}")
    print(f"Maximum y: {y_max}")
