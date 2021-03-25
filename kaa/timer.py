import time
from functools import reduce
from operator import add
from pptree import print_tree, Node

"""
Timer object containing duration data.
"""
class TimerData:

    def __init__(self, label):
        self.label = label
        self.start_time = 0
        self.end_time = 0
        self.duration = 0

    def start(self):
        self.start_time = time.time()

    def end(self):
        self.end_time = time.time()
        self.duration = self.end_time - self.start_time

"""
Static class containing all timing utility functions and statistics generating routines.
"""
class Timer:
    timer_stack = []
    timer_table = {}
    timer_child_map = {}
    timer_parent_map = {}

    @staticmethod
    def start(label):
        if not Timer.timer_stack or Timer.timer_stack[-1].label != label:
            start_timer = TimerData(label)
            start_timer.start()
            Timer.timer_stack.append(start_timer)

        else:
            raise RuntimeError("Timer of same category started previously for Timer: {}.".format(label))

    @staticmethod
    def stop(label):
        if len(Timer.timer_stack) > 0 and (current_label := Timer.timer_stack[-1].label) == label:
            end_timer = Timer.timer_stack.pop()
            end_timer.end()
            end_label = end_timer.label

            'Instantiate list for newly-found label.'
            if current_label not in Timer.timer_table:
                Timer.timer_table[end_label] = []

            Timer.timer_table[end_label].append(end_timer)

            'Populate the Timer parent/child dependencies for output later.'
            if len(Timer.timer_stack) > 0:
                if (parent_label := Timer.timer_stack[-1].label) not in Timer.timer_child_map:
                    Timer.timer_child_map[parent_label] = set()

                if current_label not in Timer.timer_child_map[parent_label]:
                    Timer.timer_child_map[parent_label].add(current_label)

                if current_label not in Timer.timer_parent_map:
                    Timer.timer_parent_map[label] = parent_label

            return end_timer.duration

        else:
            raise RuntimeError("Previous timer has not been stopped yet or timer has not been instantiated for Timer: {}.".format(label))

    @staticmethod
    def flush_timer_stack():
        Timer.timer_stack = []
        Timer.timer_table = {}
        Timer.timer_child_map = {}
        Timer.timer_parent_map = {}

    @staticmethod
    def generate_stats():
        time_data_dict = {}
        total_time = 0

        'Average time and total duration for all categories.'
        for label, times in Timer.timer_table.items():

            avg_time = Timer.avg_time(times)
            at = divmod(avg_time, 60)
            print(f"Average {label} Duration: {at[0]} min {at[1]} sec")

            tot_time = Timer.total_time(times)
            dt = divmod(tot_time, 60)
            print(f"Total {label} Duration: {dt[0]} min {dt[1]} sec \n")

            time_data_dict[label] = (at, dt, tot_time)

            if label not in Timer.timer_parent_map:
                total_time += tot_time

        'Generate Tree-like output indicating total percentage of time taken by category.'
        root = Node("Kaa Runtime")

        def construct_tree(label, parent_node):
            child_node = Node(f"{label}: {round((time_data_dict[label][2] / total_time) * 100, 2)}% \n", parent_node)
            if label not in Timer.timer_child_map:
                return

            for child_label in Timer.timer_child_map[label]:
                construct_tree(child_label, child_node)

        for label in Timer.timer_table:
            if label not in Timer.timer_parent_map:
                construct_tree(label, root)

        print_tree(root)

        Timer.flush_timer_stack()

    def avg_time(times):
        return reduce(add, [t.duration for t in times]) / len(times)

    def total_time(times):
        return reduce(add,[t.duration for t in times])
