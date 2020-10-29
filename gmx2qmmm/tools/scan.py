import numpy as np


class ReadScan:
    def __init__(self, file):
        self.file = file
        # read&match
        with open(filename, "r") as file_object:
            line = file_object.readlines()
        matchers = ["R", "A", "D"]
        matching = [s for s in line if any(xs in s for xs in matchers)]
        self.matching = np.array(matching)

        self.r_array = []
        self.a_array = []
        self.d_array = []
        for i in range(len(matching)):
            temp_line = matching[i].split()
            if (temp_line[0] == "R") and (len(temp_line) == 5):
                self.r_array = np.append(self.r_array, temp_line[1:])
            elif (temp_line[0] == "A") and (len(temp_line) == 6):
                self.a_array = np.append(self.a_array, temp_line[1:])
            elif (temp_line[0] == "D") and (len(temp_line) == 7):
                self.d_array = np.append(self.d_array, temp_line[1:])
        self.r_array = self.r_array.reshape(len(self.r_array) // 4, 4)
        self.a_array = self.a_array.reshape(len(self.a_array) // 5, 5)
        self.d_array = self.d_array.reshape(len(self.d_array) // 6, 6)

    def scanR():
        return 0

    def scan():
        return 0

    def scanR():
        return 0
