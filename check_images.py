#!/Users/karlen/anaconda2/envs/astroconda/bin/python

'''
Module for visually inspecting and interacting with UVOT images using DS9, including selection of DS9 regions for aperture photometry.
'''

import numpy as np
import os.path as path
from regions import Regions, CircleSkyRegion, CircleAnnulusSkyRegion
import warnings
import ds9_cmaps
import matplotlib.pyplot as plt

plt.style.use('dark_background')
from matplotlib import colors

from astropy.coordinates import SkyCoord
from astropy import units as u


class SourceImageViewer():
    '''Inherits from a `pyds9`_ object.
    Initialize the class with ``filepath`` input. The radii of source and background are set to suggested UVOT values.
    These are hardcoded at the moment. If they require to be changed, `aperture correction in UVOTSOURCE`_ will be necessary,
    which is not currently included in these wrappers.

    Attributes:
        filepath (str): path to image to display in DS9
        source_coords (`astropy.coordinates.SkyCoord`_): Coordinates of the source
        bkg_coords (`astropy.coordinates.SkyCoord`_): Coordinate of background region center

    .. _pyds9:
        https://github.com/ericmandel/pyds9
    .. _aperture correction in UVOTSOURCE:
        https://heasarc.nasa.gov/lheasoft/ftools/headas/uvotsource.html
    .. _astropy.coordinates.SkyCoord:
        http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord
    '''

    def __init__(self, filepath=''):
        self.filepath = filepath
        self.source_coords = None
        self.bkg_coords = None
        self.user_coords = []

    def open_fits(self):
        '''Open ``filepath``
        '''
        from astropy.io import fits
        hdulist = fits.open(self.filepath)
        self.hdu = hdulist[1]

    def draw(self):
        # self.ax.set_ylabel('Galactic Latitude')
        # self.ax.set_xlabel('Galactic Longitude')
        self.im = self.ax.imshow(self.hdu.data, origin='lower', cmap=self.cmap, norm=self.norm, zorder=1,
                                 interpolation='nearest')

        # divider = make_axes_locatable(self.ax)
        # cax = divider.append_axes('right',size="5%")
        #
        # cbar = self.fig.colorbar(self.im, cax=cax, orientation='vertical', norm=self.norm,
        #                          ticks = None, label='')

    def display_source_region(self):
        '''Display a region corresponding to source coordinates.
        '''

        self.src_handle = self.src_region.to_pixel(self.wcs).plot(ax=self.ax, zorder=2)

    def display_bkg_region(self):
        '''Tell DS9 to display a region corresponding to background center coordinates.
        '''

        self.bkg_handle = self.bkg_region.to_pixel(self.wcs).plot(ax=self.ax, zorder=2)

    def remove_regions(self, isSource):
        '''Tell DS9 to delete all regions currently open in a frame.
        '''
        if isSource:
            self.src_handle.remove()
        else:
            self.bkg_handle.remove()
        self.fig.canvas.draw()

    def load_regions(self):
        '''Load regions based on the provided ``filepath``. The `filepath`` gets parsed to
        construct the name of the region files based on the observation ID and band. If the region files
        are not found in the parent path of ``filepath``, the default coordinates will get used for
        the source (likely from SIMBAD) and for the background.
        '''

        # parse filepath to get obs ID and band
        dirpath, filename = path.split(self.filepath)
        base, extn = filename.split('_')
        obs = base[2:-3]
        band = base[-2:]

        # constructs the source and background region file paths
        regfile = path.join(dirpath, 'detect_%s_%s.reg' % (obs, band))
        bkgregfile = path.join(dirpath, 'back_%s_%s.reg' % (obs, band))

        # trying to load the background region. If the attempt fails, whatever coordinates are
        # currently stored in the object bkg_ra and bkg_dec attributes will be displayed.
        try:
            self.bkg_region = Regions.read(bkgregfile)
            self.bkg_radius = self.bkg_region.radius
            self.bkg_coords = self.bkg_region.center
        except IOError:
            print('Background region file missing. Will try default coordinates.')
            bkg_radius = 20
            self.bkg_region = CircleSkyRegion(center=self.bkg_coords, radius=bkg_radius * u.arcsec)

        self.display_bkg_region()

        # trying to load the source region. If the attempt fails, whatever coordinates are
        # currently stored in the object source_ra and source_dec attributes will be displayed.
        try:
            self.src_region = Regions.read(regfile)
            self.source_coords = self.src_region.center
            self.source_radius = self.src_region.radius
        except IOError:
            print('Region file missing. Only showing default coordinates, if any.')
            src_radius = 5
            self.src_region = CircleSkyRegion(center=self.source_coords, radius=src_radius * u.arcsec)

        self.display_source_region()

    def center_on_source(self):
        '''Center the frame on the source location.
        '''
        # centercmd = 'pan to %s %s wcs fk5' %(self.source_coords.ra.value,self.source_coords.dec.value)
        # self.set(centercmd)
        pass

    def zoom_in(self, zoom=2):
        '''Zoom in.
        '''
        self.ax.relim()
        self.ax.autoscale(True)

        center_pix = self.wcs.world_to_array_index(self.source_coords)[::-1]  # X and Y are flipped, so we flip back

        xpix_zoom = int(self.wcs.pixel_shape[1] / (zoom * 2))  # pixel_shape is flipped, [1] -> x
        ypix_zoom = int(self.wcs.pixel_shape[0] / (zoom * 2))

        self.ax.set_xlim(center_pix[0] - xpix_zoom, center_pix[0] + xpix_zoom)
        self.ax.set_ylim(center_pix[1] - ypix_zoom, center_pix[1] + ypix_zoom)

    def prettify(self):
        '''Scale the image to log and minmax and use the heat color map.
        '''
        cmap = plt.cm.get_cmap('ds9heat')
        self.cmap = cmap
        self.norm = colors.LogNorm()

    def format_frame(self):
        '''Calls ``zoom_in`` and ``prettify`` methods.
        '''
        self.zoom_in()
        self.prettify()
        self.draw()

    def update(self):
        plt.show(block=True)

    def setup_frame(self):
        '''Prepares a frame by opening a file, calling ``format_frame`` method, loading regions, and centering on the source position.
        '''

        self.open_fits()
        try:
            from astropy.wcs import WCS, FITSFixedWarning
            from astropy.visualization.wcsaxes import WCSAxes

            with warnings.catch_warnings():
                # Ignore a warning on using DATE-OBS in place of MJD-OBS
                warnings.filterwarnings('ignore', message="'datfix' made the change", category=FITSFixedWarning)
                warnings.filterwarnings('ignore', message="RADECSYS = 'FK5'", category=FITSFixedWarning)
                self.wcs = WCS(self.hdu.header)
            self.fig = plt.figure(1)
            self.ax = WCSAxes(self.fig, [0.1, 0.1, 0.8, 0.8], wcs=self.wcs)
            self.fig.add_axes(self.ax)
        except ImportError:
            self.ax = plt.subplot(111)
        self.format_frame()
        self.center_on_source()

    def find_nearest(self, array, value):
        idx = (np.abs(array - value)).argmin()
        return array[idx]

    def onclick(self, event):
        global ix, iy
        ix, iy = event.xdata, event.ydata

        # assign global variable to access outside of function
        self.user_coords.append((ix, iy))

        # Disconnect after 1 click
        if len(self.user_coords) == 1:
            self.fig.canvas.mpl_disconnect(self.cid)
            plt.close(1)
            print("Clicked")
            return

    def get_user_coords(self):
        '''Get the coordinates of the position the user clicks.
        '''

        # Simple mouse click function to store coordinates

        print("Trying user interaction")
        # Call click func
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.update()
        self.fig.canvas.draw()
        print("Clicked out")

        plt.show(block=True)

        return self.wcs.pixel_to_world(self.user_coords[-1][0], self.user_coords[-1][1])

    def prime_bkg(self):
        '''Primes the user to select a position for the center of the background region.
        An image file will be displayed and formatted. The frame will be centered on the source position and display a source region.
        DS9 will then go into interactive mode and prompt the user to select a position for the background region center on the image.
        Once the user clicks on the image, a background region will be drawn at the clicked position and the user's happiness level
        will be assessed.

        Returns:
            list: background region center RA and Dec coordinates
        '''

        self.open_fits()
        self.setup_frame()
        self.format_frame()
        self.center_on_source()

        # open loop to keep prompting the user for background selection if they're not happy after an initial attempt
        while True:
            self.setup_frame()
            self.display_source_region()
            print('Select a location for the background region.')
            self.user_coords = []
            coords = self.get_user_coords()
            bkg_ra, bkg_dec = coords.ra.deg, coords.dec.deg
            self.bkg_coords = SkyCoord('%s %s' % (bkg_ra, bkg_dec), unit=(u.deg, u.deg), frame='fk5')
            self.bkg_radius = float(input("Enter background radius in arcsec\n"))
            # self.bkg_radius = 20 #* u.arcsec
            print("User selected background ", self.bkg_coords, self.bkg_radius)
            self.update_bkg_region()
            self.setup_frame()
            self.display_source_region()
            self.display_bkg_region()

            dirpath, filename = path.split(self.filepath)
            base, extn = filename.split('_')
            obs = base[2:-3]
            band = base[-2:]

            # constructs the source and background region file paths
            bkgregfile = path.join(dirpath, 'back_%s_%s.reg' % (obs, band))

            # break
            plt.show(block=False)
            response = input('Happy with the background selection (y / anything_else)?\n')
            plt.close(1)
            if response == 'y':
                print('Great! That is all we need. Bye.')
                self.bkg_region.write(bkgregfile, format='ds9', overwrite=True)
                break
            else:
                print('Fine, try again...')
                self.remove_regions(False)

        return self.bkg_coords

    def prime_source(self):
        '''Primes the user to select a position for a source.
        An image file will be displayed and formatted. DS9 will then go into interactive mode and prompt the
        user to click on the location for the source of interest. Once the user clicks on the image, a source region
        will be drawn at the clicked position and the user's happiness level will be assessed.

        Returns:
            list: source region center RA and Dec coordinates
        '''

        self.open_fits()
        self.setup_frame()

        # open loop to keep prompting the user for source selection if they're not happy after an initial attempt
        while True:

            print('Use the plot window to select a source region.')
            coords = self.get_user_coords()
            ra, dec = coords.ra.deg, coords.dec.deg
            self.source_coords = SkyCoord('%s %s' % (ra, dec), unit=(u.deg, u.deg), frame='fk5')
            self.source_radius = float(input("Enter source radius in deg\n"))
            self.update_source_region()
            self.setup_frame()
            self.display_source_region()

            dirpath, filename = path.split(self.filepath)
            base, extn = filename.split('_')
            obs = base[2:-3]
            band = base[-2:]

            # constructs the source and background region file paths
            regfile = path.join(dirpath, 'detect_%s_%s.reg' % (obs, band))

            # break
            plt.show(block=False)
            response = input('Happy with source location selection (y/n)?\n')
            plt.close(1)
            if response == 'y':
                print('Great! That is all we need. Bye.')
                self.src_region.write(regfile, format='ds9', overwrite=True)
                break
            else:
                print('Fine, try again...')
                self.remove_regions(True)

        return self.source_coords

    def update_source_region(self):
        self.src_region = CircleSkyRegion(center=self.source_coords, radius=self.source_radius * u.arcsec,
                                          visual={'color': 'cyan'})

    def update_bkg_region(self):
        self.bkg_region = CircleSkyRegion(center=self.bkg_coords, radius=self.bkg_radius * u.arcsec,
                                          visual={'color': 'cyan'})
