# -*- coding: utf-8 -*-
# Based on the work of Lauri Kovanen, 2009-2011 (lauri.kovanen@gmail.com) 
# Department of Biomedical Engineering and Computational Science
# Aalto University School of Science

import copy

import numpy as np
import matplotlib as mpl
from matplotlib.pyplot import text
from matplotlib.patches import Rectangle
from matplotlib.colors import ColorConverter




class AlluvialDiagram:
    def __init__(self, ax, module_sizes_list, ribbon_size_matrix_list,
                 module_label_matrix=None, module_colors_list=None, 
                 ribbon_bglim=20, threshold=0.0, rainbow=False):
        """
        Plot an alluvial diagram to the ax given as parameter.
    
        Parameters
        ----------
        ax : a matplotlib.axes object
        module_sizes_list : 2D list
            a 2D list, ``module_sizes_list[i][j]'' is the size of module j in 
            time i
        ribbon_size_matrix_list : 3D list
            a 3D list, which describes how the groups change
            element with time, ``ribbon_size_matrix_list[i][j][k]``, should 
            correspond to the flow from module_sizes_list[i][j] to 
            module_sizes_list[i+1][k], each flow is a tuple, first element is 
            the outflow and second element is the inflow
        module_label_matrix : 2D list, optional
            labels for each module
        module_colors_list : iterable, optional
            colors for the first (left) modules (defaulting to gray)
        ribbon_bglim : float, optional
            if ribbon contains less than this number of nodes (or other 'mass')
            it is drawn on the background to avoid visual clutter
    
        """
        self.ax = ax
        self.ribbon_size_matrix_list = ribbon_size_matrix_list
        self.module_sizes_list = module_sizes_list
        self.module_label_matrix = module_label_matrix
        self.module_colors_list = module_colors_list
        self.ribbon_bglim = ribbon_bglim
        self.threshold = threshold
        self.rainbow = rainbow
    
    def setting_parameters(self,):
        self.ax.set_xlim([0, 2])
        self.ax.set_ylim([-0.05, 1])
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        for loc, spine in self.ax.spines.items():
            spine.set_color('none')  # don't draw spine
        if self.module_colors_list is None:
            self.module_colors_list = [["gray"] * len(module_sizes) for module_sizes in self.module_sizes_list]
    
        assert len(self.module_sizes_list) == len(self.module_colors_list)
        for sizes, colors in zip(self.module_sizes_list, self.module_colors_list):
            assert len(sizes) == len(colors)

    
        max_n_modules = np.max([len(sizes) for sizes in self.module_sizes_list])
        module_size_sum = np.max([np.sum(sizes) for sizes in self.module_sizes_list])
        self.vertical_pad_btw_modules = 0.010  # percent of
        self.vertical_pad_up_and_below = 0.00
        self.horizontal_pad_lr = 0.0
        self.individual_node_size = (1 - 2 * self.vertical_pad_up_and_below -
                                self.vertical_pad_btw_modules *
                                (max_n_modules - 1)) / module_size_sum
        self.module_width = (1 - 2 * self.horizontal_pad_lr) / (len(self.ribbon_size_matrix_list) + 1)  # should be in range [0,0.5-horizontal_pad_lr
        blank_width = (1 - 2 * self.horizontal_pad_lr) / len(self.ribbon_size_matrix_list)
        self.mw = self.module_width
        self.bw = blank_width
        # mwnsf: module width non shaded fraction
        #(to be able to differ from fully shaded to non-shaded)
        self.mwnsf = 0.1
    # plot first modules
        
    def plot_blocks(self,):
        iteration_list = [
            [module_sizes, module_colors] for module_sizes, module_colors in zip(self.module_sizes_list, self.module_colors_list)
        ]
    
        # for storing the start y coordinates
        self.module_y_starts_list = [[] for _ in  range(len(self.module_sizes_list))]
        self.module_heights_list = [[] for _ in range(len(self.module_sizes_list))]
        # plot modules
        for i, iteration_data in enumerate(iteration_list):
            module_sizes, module_colors = iteration_data
            module_y_starts = self.module_y_starts_list[i]
            module_heights = self.module_heights_list[i]
            current_y = self.vertical_pad_up_and_below
            rect_x_start = self.horizontal_pad_lr + i * (self.mw + self.bw)
            for j in range(len(module_sizes)):
                module_size = module_sizes[j]
                color = module_colors[j]
                module_y_starts.append(current_y)
                module_height = self.individual_node_size * module_size
                module_heights.append(module_height)
                rect = Rectangle((rect_x_start, current_y), self.module_width,
                             module_height, fc=color, ec="0.85")
                self.ax.add_patch(rect)
                if self.module_label_matrix is not None:
                    text(rect_x_start, current_y, self.module_label_matrix[i][j], fontsize=7)
                current_y += module_height + self.vertical_pad_btw_modules

    def plot_ribbons(self,):
        module_y_ends_list = copy.deepcopy(self.module_y_starts_list)
        curvature_param = 0.6
        # plot ribbons in order of biggest modules first?
        zorder = 0
        for t in range(len(self.ribbon_size_matrix_list)):
            for i in range(len(self.ribbon_size_matrix_list[t])):
                for j in range(len(self.ribbon_size_matrix_list[t][i])):
                    ribbon_size = self.ribbon_size_matrix_list[t][i][j]
                    if (ribbon_size[0] == 0.0) | (ribbon_size[1] == 0.0):
                        continue
                    ystart1 = self.module_y_starts_list[t][i]
                    yend1 = ystart1 + ribbon_size[0] * self.module_heights_list[t][i]
                    self.module_y_starts_list[t][i] = yend1
                    ystart2 = module_y_ends_list[t + 1][j]
                    yend2 = ystart2 + ribbon_size[1] * self.module_heights_list[t+1][j]
                    module_y_ends_list[t + 1][j] = yend2
                    # the points needed for the bezier
                    bezier_verts1 = [
                        (self.horizontal_pad_lr + t * (self.mw + self.bw) + self.module_width, ystart1),  # P0
                        (self.horizontal_pad_lr + t * (self.mw + self.bw) + self.module_width + curvature_param * self.bw, ystart1),  # P1
                        (self.horizontal_pad_lr + (t + 1) * (self.mw + self.bw) - curvature_param * self.bw, ystart2),  # P2
                        (self.horizontal_pad_lr + (t + 1) * (self.mw + self.bw), ystart2),  # P3
                        ]
                    bezier_verts2 = [
                        (self.horizontal_pad_lr + t * (self.mw + self.bw) + self.module_width, yend1),  # P0
                        (self.horizontal_pad_lr + t * (self.mw + self.bw) + self.module_width + curvature_param * self.bw, yend1),  # P1
                        (self.horizontal_pad_lr + (t + 1) * (self.mw + self.bw) - curvature_param * self.bw, yend2),  # P2
                        (self.horizontal_pad_lr + (t + 1) * (self.mw + self.bw), yend2),  # P3
                        ]
                    if max(ribbon_size) < self.ribbon_bglim:
                        use_zorder = -10000 - j
                    else:
                        use_zorder = zorder
                    if (ribbon_size[0] > self.threshold) & (ribbon_size[1] > self.threshold):
                        _plot_ribbon_using_bezier(self.ax, use_zorder, bezier_verts1,
                                          bezier_verts2,)



def _plot_ribbon_using_bezier(ax, zorder, points1, points2, color1="gray",
                              color2="gray", lw=1):
    """ Draw ribbon for alluvial diagram (see plot_alluvial)

    Parameters
    ----------
    ax : a matplotlib.axes object
    zorder : float
        the zorder for the ribbon
    points1 : iterable of float tuples
        the points, which determine the first line of the Bezier ribbon
    points2 : iterable of float tuples
        the points, which determine the second line of the Bezier ribbon
    color1 : a matplotlib compliant color definition
        color for the left side of the ribbon
    color1 : a matplotlib compliant color definition
        color for the right side of the ribbon
    lw : float
        linewidth for the bezier borders
    """
    cc = ColorConverter()
    color1 = np.array(cc.to_rgba(color1))
    color2 = np.array(cc.to_rgba(color2))
    tRange = np.linspace(0, 1, 100)
    xpointsList = []
    ypointsList = []
    for points in [points1, points2]:
        points = np.array(points)
        p1 = points[0]
        p2 = points[1]
        p3 = points[2]
        p4 = points[3]
        allPoints = (p1[:, np.newaxis] * (1 - tRange) ** 3 + p2[:, np.newaxis]
                     * (3 * (1 - tRange) ** 2 * tRange) + p3[:, np.newaxis] *
                     (3 * (1 - tRange) * tRange ** 2) + p4[:, np.newaxis] *
                     tRange ** 3)
        xpoints = allPoints[0]
        xpointsList.append(xpoints)
        ypoints = allPoints[1]
        ypointsList.append(ypoints)
        ax.plot(xpoints, ypoints, "0.85", lw=lw, zorder=zorder + 0.5)
    xpoints = xpointsList[0]
    if (mpl.colors.colorConverter.to_rgba_array(color1) ==
            mpl.colors.colorConverter.to_rgba_array(color2)).all():
        ax.fill_between(xpoints, ypointsList[0], ypointsList[1], lw=lw,
                        facecolor=color1, edgecolor=color1, zorder=zorder)
    else:
        for i in range(len(tRange) - 1):
            #mean = (tRange[i]+tRange[i+1])*0.5
            xnow = np.mean(xpoints[i:i + 2])
            norm_mean = (xnow - xpoints[0]) / (xpoints[-1] - xpoints[0])
            color = color1 * (1 - norm_mean) + color2 * norm_mean
            ax.fill_between(xpoints[i:i + 2], ypointsList[0][i:i + 2],
                            ypointsList[1][i:i + 2], lw=lw, facecolor=color,
                            edgecolor=color, zorder=zorder)                    