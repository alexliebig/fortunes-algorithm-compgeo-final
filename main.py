# CS - 6160 Computational Geometry Final Project - Fortunes Algorithm Implementation
# Creator: Alex Liebig
# Summary: In this project, my goal was to implement fortunes algorithm for creating voronoi diagrams in python only using numpy (as well other python libraries)
#          I used matplotlib to display it, as well as tkinter as a graphical user interface

import math
import numpy as np
import random
from queue import PriorityQueue
import time
from matplotlib import pyplot as plt
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

SITE = 0
CIRCLE = 1
X_DIM = 100
Y_DIM = 100
SITE_COUNT = 100

# this is the site class. This defines the x and y values of individual sites, as well as storing the voronoi vertices
class Site:
    def __init__(self, x, y, event):
        self.x = x
        self.y = y
        self.event = event
        self.voronoiVertices = []

    def __gt__(self, other):
        return other.y < self.y
    def __lt__(self, other):
        return other.y > self.y
    
    def orderVoronoiVertices(self):
        def angle(p1, p2):
            dx = p2[0] - p1[0]
            dy = p2[1] - p1[1]
            return math.atan2(dy, dx)
        try:
            bottomPoint = min(self.voronoiVertices, key=lambda p: p[1])
        except:
            return
        sorted_points = sorted(self.voronoiVertices, key=lambda p: angle(bottomPoint, p))
        self.voronoiVertices = sorted_points

# this class defines a single beach arc. The x intercepts are stored here, and the site it is referencing is dependency injected into it
class BeachArc:
    def __init__(self, site):
        self.site = site
        self.left_intercept_x = -np.inf
        self.right_intercept_x = np.inf
        self.left_point_queue = []
        self.right_point_queue = []
        self.vertex_queue = []

# directrix is the sweep line
# focus is the current site
# x is the site point
# formula: a is focus x, b is focus y
#    ((x-a)^2+b^2-c^2)
# y = ----------------
#         2(b-c)
def calculate_y_parabola(x, focus, directrix):
   numerator = (x-focus[0])**2 + focus[1]**2 - directrix**2
   denominator = 2 * (focus[1] - directrix)
   if denominator == 0:
      return x
   return numerator / denominator

# used to calculate x intercept of parabolas
# sets both focus parabolas equal to eachother
def calculate_x_parabola_intercept(focus1, focus2, directrix):
    D0 = (focus1[0]**2 + focus1[1]**2 - directrix**2)
    D1 = (focus2[0]**2 + focus2[1]**2 - directrix**2)
    E0 = 4 * (focus1[1] - directrix)
    E1 = 4 * (focus2[1] - directrix)

    if E0 == 0:
        return focus1[0]
    if E1 == 0:
        return focus2[0]

    a1 = 1 / E0
    b1 = -2 * focus1[0] / E0
    c1 = D0 / E0
    
    a2 = 1 / E1
    b2 = -2 * focus2[0] / E1
    c2 = D1 / E1
    
    # get a b and c values for quadratic formula
    a = a1 - a2
    b = b1 - b2
    c = c1 - c2

    if b**2 - (4*a*c) < 0:
        return -1

    x1 = (-b - math.sqrt(b**2 - (4*a*c)))/(2*a)
    x2 = (-b + math.sqrt(b**2 - (4*a*c)))/(2*a)
    if x1 >= focus1[0] and x1 <= focus2[0]:
       return x1
    return x2

# this simply just creates random points to plot on the voronoi diagram
def createRandomSites(n):
    sites = [None]*n
    for i in range(n):
        sites[i] = Site(random.random() * X_DIM, random.random() * Y_DIM, SITE)

    # adds extra points to make voronoi diagram look decent
    sites.append(Site( X_DIM*2, Y_DIM*2*1.00, SITE))
    sites.append(Site( X_DIM*2,-Y_DIM*2*1.01, SITE))
    sites.append(Site(-X_DIM*2, Y_DIM*2*1.10, SITE))
    sites.append(Site(-X_DIM*2,-Y_DIM*2*1.11, SITE))
    return sites

# insert into list based on x value
def insert_into_active_list(active_set, curr_point, directrix):
    if len(active_set) == 0:
        active_set.append(BeachArc(curr_point))
        return 0
    else: 
        curr_index = 0
        if len(active_set) == 1:
            active_set.insert(curr_index + 1, BeachArc(active_set[curr_index].site))
            active_set.insert(curr_index + 1, BeachArc(curr_point))
            updateIntersections(active_set, directrix)
            return 1
        else:
            updateIntersections(active_set, directrix)
            for i,val in enumerate(active_set):
                if curr_point.x > val.left_intercept_x and curr_point.x < val.right_intercept_x:
                    curr_index = i
                    break

            # break everything up
            originalPoint = active_set[curr_index].site
            active_set.insert(curr_index, BeachArc(curr_point))
            active_set.insert(curr_index, BeachArc(originalPoint))
            return curr_index+1

# this function updates each intersections of the beach arc in the beach line
def updateIntersections(active_set, directrix):
    for i,val in enumerate(active_set):
        if i == 0:
            continue
        intercept = calculate_x_parabola_intercept([active_set[i-1].site.x, active_set[i-1].site.y], [active_set[i].site.x, active_set[i].site.y], directrix)
        active_set[i-1].right_intercept_x = intercept
        active_set[i].left_intercept_x = intercept

# this function gets the bottom of the circle
def circle_bottom(x1, y1, x2, y2, x3, y3):
    a = x1 - x2
    b = y1 - y2
    c = x1 - x3
    d = y1 - y3
    e = (x1**2 - x2**2) + (y1**2 - y2**2)
    f = (x1**2 - x3**2) + (y1**2 - y3**2)
    den = 2 * (a*d - b*c)
    if den == 0:
        return False
    cx = (d*e - b*f) / den
    cy = (a*f - c*e) / den
    r = math.sqrt((x1 - cx)**2 + (y1 - cy)**2)
    return cy - r

#this function processes the circle events that are encountered
def processCircleEvent(i, Q, active_set, directrix):
    leftArc = active_set[i-1] if i-1 >= 0 else None
    rightArc = active_set[i+1] if i+1 < len(active_set) else None
    if leftArc is None:
        return None
    if rightArc is None:
        return None
    if hash(leftArc.site) == hash(rightArc.site):
        return None
    
    leftSite = leftArc.site
    rightSite = rightArc.site
    currSite = active_set[i].site

    bottomVal = circle_bottom(leftSite.x, leftSite.y, currSite.x, currSite.y, rightSite.x, rightSite.y)
    if bottomVal < directrix:
        newPoint = Site(None, bottomVal, CIRCLE)
        Q.put(((Y_DIM - bottomVal), newPoint))

# this is where the MAGIC happens ;). This is where the voronoi diagram is computed!
def plot_voronoi(points):
    # create priority queue
    Q = PriorityQueue()
    active_set = []
    for point in points:
        Q.put(((Y_DIM - point.y), point))

    while(Q.empty() != True):
        point = Q.get()[1]
        directrix = point.y
        if point.event == SITE:
            # insert into list
            i = insert_into_active_list(active_set, point, directrix)
            processCircleEvent(i-1, Q, active_set, directrix)
            processCircleEvent(i+1, Q, active_set, directrix)
        elif point.event == CIRCLE:
            # adjust all intersections in active list
            updateIntersections(active_set, directrix)
            # remove site if left and right intersection points are equal
            for i,beachArc in enumerate(active_set):
                zeroVal = math.fabs(beachArc.left_intercept_x - beachArc.right_intercept_x)
                if math.isclose(0, zeroVal, rel_tol=1e-6, abs_tol=1e-6):
                    voronoiVertex = [beachArc.left_intercept_x, calculate_y_parabola(beachArc.left_intercept_x, [beachArc.site.x,beachArc.site.y], directrix)]
                    active_set[i-1].site.voronoiVertices.append(voronoiVertex)
                    active_set[ i ].site.voronoiVertices.append(voronoiVertex)
                    active_set[i+1].site.voronoiVertices.append(voronoiVertex)
                    active_set.pop(i)
                    processCircleEvent(i-1, Q, active_set, directrix)
                    processCircleEvent(i, Q, active_set, directrix)
                    break

# plots the voronoi polygons
def plotVoronoiPolygons(points):
    for point in points:
        point.orderVoronoiVertices()
        xVals = np.array(point.voronoiVertices)[:,0]
        yVals = np.array(point.voronoiVertices)[:,1]
        if len(xVals) >= 3:
            plt.fill(xVals, yVals, alpha=0.5)

# this is the main event loop for the voronoi diagram that implements fortunes algorithm
def voronoiEventLoop():
    SITE_COUNT = int(entry.get())
    while(1):
        points = createRandomSites(SITE_COUNT)
        ax.clear()
        plt.title("Fortune's Algorithm Voronoi Diagram Viewer")
        for point in points:
            plt.plot(point.x,point.y,'bo', markersize=4)
        try:
            plot_voronoi(points)
        except:
            ax.clear()
            continue
        plt.xlim([0, X_DIM])
        plt.ylim([0, Y_DIM])
        plotVoronoiPolygons(points)
        break
    canvas.draw()

# This functions calls the voronoi event loop from the Enter keyword
def on_enter_key(event):
    voronoiEventLoop()

# this function runs benchmarks using the "T" key from the window
def on_benchmark_key(event):
    siteCounts = [50, 100, 200, 500, 1000, 2000, 5000]
    times = np.zeros(7)
    for i,siteCount in enumerate(siteCounts):
        points = createRandomSites(siteCount)
        startTime = None
        endTime = None
        while(1):
            try:
                startTime = time.time()
                plot_voronoi(points)
                endTime = time.time()
            except:
                continue
            break
        print(f"Site Count {siteCount} total time in seconds {endTime - startTime}")
        times[i] = endTime - startTime
    
    # Save Benchmarks as a plot
    secondaryFig = Figure(figsize=(5, 4), dpi=500)
    secondaryAx = secondaryFig.add_subplot(111)
    secondaryAx.plot(siteCounts, times, color='red')
    secondaryAx.set_xlabel("Site Count")
    secondaryAx.set_ylabel("Time to Execute (in seconds)")
    secondaryAx.set_title("Performance of Fortune's Algorithm")
    secondaryFig.savefig("performance.png")

if __name__ == "__main__":
    # Setup Tkinter Interface
    root = tk.Tk()
    root.title("Fortune's Algorithm Voronoi Plot Viewer")
    root.geometry("500x500")

    # initial set-up
    while(1): # used in case there's an exception
        points = createRandomSites(SITE_COUNT)
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_aspect('equal')
        plt.title("Fortune's Algorithm Voronoi Diagram Viewer")
        for point in points:
            plt.plot(point.x,point.y,'bo', markersize=4)
        try: # wayt to implement exception handling
            plot_voronoi(points)
        except:
            ax.clear()
            continue
        plt.xlim([0, X_DIM])
        plt.ylim([0, Y_DIM])
        plotVoronoiPolygons(points)
        break

    # Create Canvas
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    # Set-up input field
    validate_int = root.register(lambda s: s.isdigit())
    entry = tk.Entry(root, validate="key", validatecommand=(validate_int, '%S'), width=int(100*0.5), textvariable=100)
    entry.insert(0, '100')
    entry.pack(side=tk.LEFT)

    # Create a submit button to update the plot
    submit = tk.Button(root, text="Randomize Points", command=voronoiEventLoop, width=int(100*0.5))
    submit.pack(side=tk.RIGHT)

    # Handle Closing Event
    root.protocol("WM_DELETE_WINDOW", lambda: (root.destroy(), root.quit()))

    # Handle Enter Event
    root.bind("<Return>", on_enter_key)
    root.bind("<t>", on_benchmark_key)

    print("Usage: type a valid integer amount for the number of sites you want to display.")
    print("Usage: You can hit <Enter> or click \"Randomize Points\" to display the new voronoi diagram.")
    print("Usage: You can press <T> to run benchmarks on the algorithm.")

    # Tkinter Event Loop Start
    tk.mainloop()