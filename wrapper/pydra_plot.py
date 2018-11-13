#!/usr/bin/env python3
#
# JM, 30 Oct 2018
#
# some plotting functions relating to pydra

from pydra_analysis import *
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["mathtext.rm"] = "serif"
plt.rcParams["image.cmap"] = "RdBu_r" # "*_r" is reverse of standard colour
plt.rcParams["axes.formatter.limits"] = [-3, 3]

#-------------------------------------------------------------------------------
# plot gamma_m and phi_m

def plot_geom_param(data_geom, y_vec, t_now):

  # momentum things
  ax = plt.subplot2grid((4, 4), (0, 0), colspan = 1)
  ax.plot(data_geom.gamma_m_L1L2[:, 1], y_vec)
  ax.set_xlim([0, 1.1])
  ax.set_title(r"$\gamma_{m, 1} (t = %g)$" % t_now)
  ax.grid()

  ax = plt.subplot2grid((4, 4), (1, 0), colspan = 1)
  ax.plot(data_geom.gamma_m_L1L2[:, 0], y_vec)
  ax.set_xlim([0, 1.1])
  ax.set_title(r"$\gamma_{m, 2} (t = %g)$" % t_now)
  ax.grid()

  ax = plt.subplot2grid((4, 4), (0, 2), colspan = 1)
  ax.plot(data_geom.phi_m_L1L2[:, 1] / np.pi, y_vec)
  ax.set_xlim([-0.6, 0.6])
  ax.set_title(r"$\phi_{m, 1} / \pi (t = %g)$" % t_now)
  ax.grid()

  ax = plt.subplot2grid((4, 4), (1, 2), colspan = 1)
  ax.plot(data_geom.phi_m_L1L2[:, 0] / np.pi, y_vec)
  ax.set_xlim([-0.6, 0.6])
  ax.set_title(r"$\phi_{m, 2} / \pi (t = %g)$" % t_now)
  ax.grid()

  ax = plt.subplot2grid((4, 4), (0, 1), colspan = 1)
  ax.plot(data_geom.gamma_m_btbc[:, 0], y_vec)
  ax.set_xlim([0, 1.1])
  ax.set_yticklabels([])
  ax.set_title(r"$\gamma_{m, \rm{bt}} (t = %g)$" % t_now)
  ax.grid()

  ax = plt.subplot2grid((4, 4), (1, 1), colspan = 1)
  ax.plot(data_geom.gamma_m_btbc[:, 1], y_vec)
  ax.set_xlim([0, 1.1])
  ax.set_yticklabels([])
  ax.set_title(r"$\gamma_{m, \rm{bc}} (t = %g)$" % t_now)
  ax.grid()

  ax = plt.subplot2grid((4, 4), (0, 3), colspan = 1)
  ax.plot(data_geom.phi_m_btbc[:, 0] / np.pi, y_vec)
  ax.set_xlim([-0.6, 0.6])
  ax.set_title(r"$\phi_{m, \rm{bt}} / \pi (t = %g)$" % t_now)
  ax.grid()

  ax = plt.subplot2grid((4, 4), (1, 3), colspan = 1)
  ax.plot(data_geom.phi_m_btbc[:, 1] / np.pi, y_vec)
  ax.set_xlim([-0.6, 0.6])
  ax.set_title(r"$\phi_{m, \rm{bc}} / \pi (t = %g)$" % t_now)
  ax.grid()
  
  # buoyancy things
  ax = plt.subplot2grid((4, 4), (2, 0), colspan = 1)
  ax.plot(data_geom.gamma_b[:, 1], y_vec)
  ax.set_title(r"$\gamma_{b, 1} (t = %g)$" % t_now)
  ax.grid()

  ax = plt.subplot2grid((4, 4), (3, 0), colspan = 1)
  ax.plot(data_geom.gamma_b[:, 0], y_vec)
  ax.set_title(r"$\gamma_{b, 2} (t = %g)$" % t_now)
  ax.grid()

  ax = plt.subplot2grid((4, 4), (2, 1), colspan = 1)
  ax.plot(data_geom.phi_b / np.pi, y_vec)
  ax.set_xlim([-1.1, 1.1])
  ax.set_title(r"$\phi_{b} / \pi (t = %g)$" % t_now)
  ax.grid()

  ax = plt.subplot2grid((4, 4), (2, 2), colspan = 1)
  ax.plot(data_geom.gamma_t[:, 1], y_vec)
  ax.set_xlim([0, 1.1])
  ax.set_title(r"$\gamma_{t, 1} (t = %g)$" % t_now)
  ax.grid()

  ax = plt.subplot2grid((4, 4), (3, 2), colspan = 1)
  ax.plot(data_geom.gamma_t[:, 0], y_vec)
  ax.set_xlim([0, 1.1])
  ax.set_title(r"$\gamma_{t, 2} (t = %g)$" % t_now)
  ax.grid()

  ax = plt.subplot2grid((4, 4), (2, 3), colspan = 1)
  ax.plot(data_geom.phi_t[:, 1] / np.pi, y_vec)
  ax.set_xlim([-0.1, 0.6])
  ax.set_title(r"$\phi_{t, 1} / \pi (t = %g)$" % t_now)
  ax.grid()

  ax = plt.subplot2grid((4, 4), (3, 3), colspan = 1)
  ax.plot(data_geom.phi_t[:, 0] / np.pi, y_vec)
  ax.set_xlim([-0.1, 0.6])
  ax.set_title(r"$\phi_{t, 2} / \pi (t = %g)$" % t_now)
  ax.grid()
  
#TODO: write a subfunction to plot the processed (e.g. alpha, sin(phi_b)) data  
