#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created by Giorgia Bussu on february 16 2022 in code form
    adapting the original pop-out visual attention task after sound cue created by Lowe Wilsson
    using PsychoPy3 Experiment Builder (v2021.2.3),
    on februari 08, 2022, at 15:31
If you publish work using this script the most relevant publication is:

    Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
        PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
        https://doi.org/10.3758/s13428-018-01193-y

"""

from __future__ import absolute_import, division

from psychopy import locale_setup
from psychopy import prefs
prefs.hardware['audioLib'] = 'ptb'
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)

import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle, choice as randchoice
import os  # handy system and path functions
import sys  # to get file system encoding
import pandas as pd

from psychopy.hardware import keyboard

### BEGIN SET CONSTANTS ###
# unless otherwise specified,
# stimuli sizes/coordinates are specified
# in degrees (deg), and times are specified
# in seconds

# notes about 'Kleberg 2019' refer to
# "How Infants’ Arousal Influences Their Visual Search"
# - Kleberg, del Bianco, Falck-Ytter 2019

# variable for activating (True) or disabling (False)
# debug mode. make sure that this is set to 
# **False**, or the variable and associated debugging
# code is entirely removed, before running the experiment
# with real participants.
# note that for the time being, **debug mode will crash**
# after 4 trials unless you also reduce the trial_loop to just
# 4 iterations. the reason is that less images are loaded, to
# decrease startup time.
DEBUG_ON = False

# eyetracker update frequency (unfortunately this has to be updated
# both here and in YAML file), eg 600 if eyetracker captures
# gaze 600 times/second (600Hz).
EYETRACKER_HZ = 600

# load screen image width/height (this is specified
# separately from other images, since it doesn't
# share specifications with anything)
LOAD_IMG_WIDTH = 15.65
LOAD_IMG_HEIGHT = 12


# visual stimuli distance from centre of screen,
# in x/y direction
IMG_Y_OFFSET = 6
IMG_X_OFFSET = 12

# 'in-trial' visual stimuli max width/heights
# (stimuli are scaled, maintaining proportions,
# to being as large as possible without surpassing
# the maximum width/height)
IMG_MAX_WIDTH = 9
IMG_MAX_HEIGHT = 9
# attention grabber max width/heights
ATT_MAX_WIDTH = 5
ATT_MAX_HEIGHT = 5

# each visual stimulus also has an associated area of interest
# (AOI; see trial code components for details on implementations
# of this), with its own (same as or larger than the stimulus) 
# width/height.
# trial AOI
TRIAL_AOI_WIDTH = 10
TRIAL_AOI_HEIGHT = 10
# attention grabber AOI
ATT_AOI_WIDTH = 10
ATT_AOI_HEIGHT = 10

# trial audio low/high volume levels. note that these describe
# *PsychoPy* sound levels, which are directly proportional to the
# physical sound/air pressure. this is different from decibels -
# look for info online if you need to understand the relationship.
TRIAL_AUDIO_1 = 0
TRIAL_AUDIO_2 = 0.03125 #0.286
TRIAL_AUDIO_3 = 0.0625 #0.429
TRIAL_AUDIO_4 = 0.125 #0.571
TRIAL_AUDIO_5 = 0.25 #0.714
TRIAL_AUDIO_6 = 0.5 #0.857
TRIAL_AUDIO_7 = 1

# path to csv file which specifies image file paths,
# original image widths/heights (in pixels), type
# of image (social/geometric shape/man-made/natural/attention grabber),
# and sex (male/female)
IMG_SPEC_PATH = (
    'stimuli_specifications/visual_stimuli_specifications.csv'
)
# if debug mode is on, use less images to decrease load time
if DEBUG_ON:
    IMG_SPEC_PATH = (
        'stimuli_specifications/visual_stimuli_specifications_DEBUG.csv'
    )

# path to csv file which specifies sound file paths,
# sound types (social/nonsocial) and sex (male/female/
# nonsocial)
AUDIO_SPEC_PATH = (
    'stimuli_specifications/auditory_6stimuli_specifications.csv'
)
# path to csv file which specifies attention grabber
# image file paths and original image 
# widths/heights (in pixels)
ATT_IMG_SPEC_PATH = (
    'stimuli_specifications/'
    'attentiongrabber_images_specifications.csv'
)
# path to csv file which specifies 
# attention grabbing sound file paths
ATT_AUDIO_SPEC_PATH = (
    'stimuli_specifications/'
    'attentiongrabber_audio_specifications.csv'
)

# attention grabber earliest/latest time of
# onset (counting from trial start time) 
# -- onset directly at trial start, based on 
# internal instruction from Falck-Ytter and Da Silva-- 
ATT_GRAB_EARLIEST = 0
ATT_GRAB_LATEST = 0

# attention grabber __maximum__ duration (ie if eyetracker is malfunctioning
# and/or the participant is not at all looking at the attention grabber,
# force experiment to continue after this duration)
ATT_GRAB_MAX_DUR = 5

# sound min/max delay after gaze is directed
# at attention grabber
# -- onset directly after gaze capture, based on 
# internal instruction from Falck-Ytter and Da Silva-- 
AUDITORY_MIN_DELAY = 0
AUDITORY_MAX_DELAY = 0

AUDITORY_MIN_DELAY_HABIT = 1
AUDITORY_MAX_DELAY_HABIT = 1

# min/max delay from sound to
# visual stimuli onset 
# -- based on Kleberg 2019 --
VISUAL_MIN_DELAY = 0.08
VISUAL_MAX_DELAY = 0.4
# visual stimuli 
# min/max presentation duration
# -- based on Kleberg 2019 --
VISUAL_MIN_DUR = 3
VISUAL_MAX_DUR = 3

# min/max duration of 'blank' screen after
# visual stimuli have finished (and before
# trial end)
# -- based on internal instruction from Falck-Ytter 
# and Da Silva--
BLANK_MIN_DUR = 0.5
BLANK_MAX_DUR = 0.5

# how much time (in seconds), after
# experimenter has played attention grabbing sound,
# must pass before experimenter is allowed to play
# a new attention grabbing sound
ATT_SOUND_COOLDOWN_TIME = 2

## habituation video specifications ##

# habituation video min/max delay after gaze is directed
# at attention grabber
# -- onset directly after gaze capture, based on 
# internal instruction from Falck-Ytter and Da Silva-- 
VIDEO_MIN_DELAY = 0
VIDEO_MAX_DELAY = 0

VIDEO_MIN_DUR = 4
VIDEO_MAX_DUR = 4
# ATTENTION GRABBER AND BLANK END SCREEN SPECIFICS ARE THE SAME AS ABOVE (EXPERIMENT 1, POPOUT)

# path to csv file which specifies habituation
# video file paths and color/type
# randomized each infant in terms of green/red sync/nosync

carsequence_rand = np.random.choice(np.arange(100))  

# A sequence: red sync, green no sync
if carsequence_rand % 2 == 0: 
    
    HVIDEO_SPEC_PATH = ('stimuli_specifications/habituation_Astimuli_specifications.csv')
    HAUDIO_SPEC_PATH = ('stimuli_specifications/habituation_sounds_Astimuli_specifications.csv')

# B sequence: red no sync, green sync
if carsequence_rand % 2 == 1: 
     HVIDEO_SPEC_PATH = ('stimuli_specifications/habituation_Bstimuli_specifications.csv')
     HAUDIO_SPEC_PATH = ('stimuli_specifications/habituation_sounds_Bstimuli_specifications.csv')

# generating random sequence for order sync/noync presentation of habituation videos
habit_trial_list = np.arange(6)
shuffle(habit_trial_list)

a =  0

HVIDEO_WIDTH_DEG=28
# the width/height aspect ratio is the same for all pause videos
# being used, 1920(width):1080(height).
HVIDEO_HEIGHT_DEG=1080/1440 * HVIDEO_WIDTH_DEG

## end pause video specifications ##

## attention grabber manipulations ##
# --- pulsating attention grabber:
# period, ie time spent on a whole
# shrink/grow cycle
PULSE_PERIOD = 1
# relative amplitude - for example, 0.2
# means that the attention grabber will
# grow to 120% and shrink to 80%,
# 0.3 -> 130%/70%
PULSE_AMP = 0.4
# --- flickering attention grabber:
# period, ie time spent on a whole
# fade-out/fade-in cycle
FLICKER_PERIOD = 0.7
# least opacity - for example, 0.5
# means that the attention grabber will
# fade out to 50% opacity/'visibility',
# then fade back in again to 100%
FLICKER_LEAST_OPACITY = 0.5
# --- bouncing attention grabber:
# period, ie time spent on a whole
# 'jump', going up/down and back to
# original position
BOUNCE_PERIOD = 1
# height, ie how many degrees (unless
# experiment unit is changed to something
# else than degree) the attention grabber
# should move up/down
BOUNCE_HEIGHT = 0.5

## end attention grabber manipulations ##

## messages to participant ##
# load screen message.
# "Please wait for a bit..."
LOAD_MESSAGE = (
    "Vänligen vänta lite..."
)
# end screen message.
# "With that, the experiment is finished."
# "Thank you!"
END_MESSAGE = (
    "Nu är experimentet klart.\n"
    "Tack!"
)

# text height (in degrees) for messages to participant
MESSAGE_TXT_HEIGHT = 1

## end messages to participant ##

### END SET CONSTANTS ###
# load necessary libraries
from psychopy.iohub import launchHubServer

# set constant that specifies if the eyetracker
# is to be 'mocked', ie the mouse is used instead
# of an actual eyetracker
MOCK_ON = False

# import the pandas package (which is already being 
# loaded and used by PsychoPy in the background, 
# so there's no performance loss)
# in order to be able to read and process
# CSV files with stimuli specifications
import pandas as pd
# read in image/audio stimuli specifications
# and put them in pandas dataframes
# 'in-trial' visual stimuli
img_df = pd.read_csv(IMG_SPEC_PATH)
# 'in-trial' audio stimuli
audio_df = pd.read_csv(AUDIO_SPEC_PATH)
# attention grabbing audio stimuli
att_audio_df = pd.read_csv(ATT_AUDIO_SPEC_PATH)
# attention grabber image specifications
att_img_df = pd.read_csv(ATT_IMG_SPEC_PATH)

def get_trial_img_ls():
    """
    Randomly grabs
    1 social image stimulus component
    1 geometrical shape stimulus component
    1 natural object stimulus component
    1 man-made object stimulus component
    and returns them as a randomly-ordered
    4-element numpy array.
    Note that this function relies on the lists `img_ls_social`
    and `img_ls_nonsocial` already having been created before
    this function is called. See the 'Begin Experiment' tab,
    where these two lists are created.
    """
    # 'pop off' an image from each category's list of images. this means
    # that, for each category, the last image is 'plucked' from the end of 
    # the list of images, and 'given to' the image variable to the left
    # of the '=' operator/sign. doing this ensures that once an image has
    # been fetched for use in one trial, it won't be used for another
    soc_img = img_ls_social.pop()
    geo_img = img_ls_geo.pop()
    manm_img = img_ls_manm.pop()
    nat_img = img_ls_nat.pop()
    # combine the four images into a list
    trial_img_ls = [soc_img, geo_img, manm_img, nat_img]
    return trial_img_ls

def calc_width_height(orig_w, orig_h, max_w, max_h):
    """
    Calculates PsychoPy image stimulus width/height specifications based on
    original image's width/height (orig_w, orig_h) and the maximum width/height
    in PsychoPy.
    """
    w_ratio = max_w / orig_w
    h_ratio = max_h / orig_h
    if (w_ratio < h_ratio):
        return w_ratio * orig_w, w_ratio * orig_h
    else:
        return h_ratio * orig_w, h_ratio * orig_h

def random_20thsec(max_n_20ths):
    """
    Randomly picks a number in the range [0, <max_n_20ths>*0.05],
    in order to randomly pick a duration using increments of 0.05s.
    Note that the range is inclusive on both ends, ie if you pass
    in 4 you might get 0, 0.05, 0.10, 0.15 or 0.2 back.
    """
    num_increments_0p05 = np.random.randint(max_n_20ths+1)
    return num_increments_0p05/20

def generate_img_components_from_df(
    df, 
    max_width=IMG_MAX_WIDTH, 
    max_height=IMG_MAX_HEIGHT,
    img_category="default",
    with_lin_interpolation=True
):
    """
    (used in Begin Experiment) Takes in a data frame (df), describing images, with the following columns:
        - 'width', specifying image's original width
        - 'height', specifying image's original height
        - 'file_path', specifying image's file path (relative to .psyexp file, or absolute path)
    Returns a list of PsychoPy ImageStim instances, ie image components, which are
    based on the passed image widths/heights/paths, except down/up-sized based on
    the experiment's maximum allowed image width/heights (IMG_MAX_WIDTH/IMG_MAX_HEIGHT).
    
    Optionally, a (string) img_category argument can be passed to this 
    function (defaults to "default"), which will be set as the 'img_category'
    property of all resulting image components. Eg if 'manmade' is passed as
    an argument, then all produced image components will have an img_category 
    property with value 'manmade'.
    
    'with_lin_interpolation' is used to specify if linear interpolation should be
    used. This is necessary for photo-based image components (see 'Begin Experiment' 
    for more information on a bug that has to be handled, however) but should not be used
    for simpler images, such as the attention grabbers which are non-photo-based.
    """
    img_ls = []
    # loop over all rows of the image dataframe,
    # creating an image component with each
    # row's specifications, and appending it to
    # the list
    for row_ind in range(len(df)):
        row = df.iloc[row_ind]
        width_height = calc_width_height(
            row['width'],
            row['height'],
            max_width,
            max_height
        )
        new_img = visual.ImageStim(
            win=win,
            name=f'image_social_{row_ind}', 
            image=df.iloc[row_ind]['file_path'], 
            mask=None,
            ori=0, 
            pos=(0, 0), 
            size=width_height,
            color=[1,1,1], 
            colorSpace='rgb', 
            opacity=1,
            flipHoriz=False, 
            flipVert=False,
            texRes=2048, 
            interpolate=with_lin_interpolation, 
            depth=0.0
        )
        new_img.color='white'
        if img_category:
            new_img.img_category = img_category
        img_ls.append(new_img)
    # shuffle (randomly order) 
    # the generated list of image components,
    # to make sure that each trial run gets a
    # completely random order of images
    shuffle(img_ls)
    # return the generated list of image components
    return img_ls

def generate_position_indices(previous_indices = []):
    """
    Generates visual stimuli position indices. eg 0 might correspond to
    top left, 1 top right, etc. The interpretation of the indices
    depends on how the generated values are used, since this
    function doesn't know about the actual positions, it just generates
    numbers. See the 'Begin Routine' tab for code where this function is used.
    
    As the first argument, 'previous_positions', this function takes in
    a list of position indices used in the previous (if any) trial.
    The function ensures that the new indices it returns doesn't match
    the previous trials' at any point in the list. So for example,
    if this function is passed [0, 2, 3, 1], then it could return
    [1, 3, 0, 2], but not [1, 4, 2, 3], since the 1 would be at the start
    of the list in both 'previous_indices' and the in the returned list.
    """
    # if the passed list of previous indices is empty (ie there hasn't been a
    # previous trial), fill it up with numbers that won't cause a collision
    if not previous_indices:
        previous_indices = [-1, -1, -1, -1]
    # get the length of the 'previous_indices' list, to make sure
    # that as many indices are generated here
    num_inds = len(previous_indices)
    # generate a list of the possible index numbers
    possible_indices = list(range(num_inds))
    
    found_valid_list = False
    # keep generating lists of position indices until one that doesn't
    # 'collide' with the passed previous_indices list is found
    while not found_valid_list:
        collides = False
        # randomly reorder the list of possible index numbers
        shuffle(possible_indices)
        # check for any collisions
        for i in range(num_inds):
            if possible_indices[i] == previous_indices[i]:
                collides = True
        if not collides:
            found_valid_list = True
    return possible_indices

def assign_volume(sound_component):
    """
    Takes in a sound component and increments volume level by 1 until max level. 
    """
    sound_option_ls = [TRIAL_AUDIO_1, TRIAL_AUDIO_2, TRIAL_AUDIO_3, TRIAL_AUDIO_4, TRIAL_AUDIO_5, TRIAL_AUDIO_6, TRIAL_AUDIO_7]
    if sound_component.volume == TRIAL_AUDIO_1:
        sound_component.volume = TRIAL_AUDIO_2
        
    elif sound_component.volume == TRIAL_AUDIO_2:
        sound_component.volume = TRIAL_AUDIO_3
        
    elif sound_component.volume == TRIAL_AUDIO_3:
        sound_component.volume = TRIAL_AUDIO_4
        
    elif sound_component.volume == TRIAL_AUDIO_4:
        sound_component.volume = TRIAL_AUDIO_5
        
    elif sound_component.volume == TRIAL_AUDIO_5:
        sound_component.volume = TRIAL_AUDIO_6
        
    elif sound_component.volume == TRIAL_AUDIO_6:
        sound_component.volume = TRIAL_AUDIO_7
        

# this code snippet defines functions which are
# to be used with attention grabber image components,
# (in separate code snippets' 'each frame' tab)
# to animate them

def pulsate_image(
    img, 
    time_val, 
    orig_size,
    rel_amplitude=PULSE_AMP, 
    period=PULSE_PERIOD
):
    """
    Resizes a passed image component based on a passed
    time (in seconds) and a passed relative amplitude.
    E. g. if the passed relative amplitude is 0.2,
    if this function is called continuously with the
    image component then the image will pulsate
    between 80% and 120% of its original size.
    :param img: ImageStim - image component to be resized/pulsate
    :param time_val: float - current time in seconds (ie the 't'  variable in PsychoPy Builder)
    :param orig_size: 2-element tuple - describes the image's original size, eg (0.3, 0.4)
    :param rel_amplitude: float - relative amplitude, see above
        (defaults to value set in 'setup' routine's code component)
    :param period: float - seconds per complete pulse/cycle
        (defaults to value set in 'setup' routine's code component)
    """
    # using the numpy package `sin`, sine, function here. the sine function is ideal for
    # describing cyclical processes. for information about the function and the
    # mathematics involved, see eg https://arrayjson.com/numpy-sin/
    # also, the numpy 'pi' variable, which holds an 
    # approximation of the mathematical constant pi, is used
    resize_factor = rel_amplitude * sin(time_val / period * 2 * pi) + 1
    new_size = tuple((resize_factor * x for x in orig_size))
    img.size = new_size

def flicker_image(
    img, 
    time_val, 
    orig_opacity=1,
    least_opacity=FLICKER_LEAST_OPACITY, 
    period=FLICKER_PERIOD
):
    """
    Makes a passed image component 'flicker', ie
    fade in and out, based on a passed
    time (in seconds) and a passed least opacity
    value.
    E. g. if the passed least opacity value is 0.6,
    if this function is called continuously with the
    image component then the image will pulsate
    between 60% and 100% opacity.
    :param img: ImageStim - image component to make fade in/out
    :param time_val: float - current time in seconds (ie the 't'  variable in PsychoPy Builder)
    :param orig_opacity: float - describes the image's original opacity (defaults to 1, ie 100% visibility)
    :param least_opacity: float - least opacity, see above 
        (defaults to value set in 'setup' routine's code component)
    :param period: float - seconds per complete fade-in/fade-out cycle 
        (defaults to value set in 'setup' routine's code component)
    """
    # find 'middle' opacity, eg 0.75 if original opacity is 1 and least_opacity is
    # set to 0.5
    middle_opacity = (orig_opacity + least_opacity)/2
    # find the difference between 'middle' and original opacity, eg 0.25 in the example above
    mid_to_orig_diff = orig_opacity - middle_opacity
    # for info on `sin`, see comment in `pulsate_image` code above
    new_opacity = middle_opacity + mid_to_orig_diff * sin(time_val / period * 2 * pi)
    img.opacity = new_opacity

def bounce_image(
    img, 
    time_val, 
    orig_y_coord=0,
    bounce_height=BOUNCE_HEIGHT, 
    period=BOUNCE_PERIOD
):
    """
    Makes a passed image component 'bounce', ie
    move up and down, based on a passed
    time (in seconds) and a passed bounce height
    value.
    E. g. if the passed bounce height value is 1,
    if this function is called continuously with the
    image component then the image will continuously
    move up 1 degree, go back to original position, go down 1 deg,
    go back to original position, and so on.
    :param img: ImageStim - image component to make bounce
    :param time_val: float - current time in seconds (ie the 't'  variable in PsychoPy Builder)
    :param orig_y_coord: float - describes the image's original y coordinate/position 
        (defaults to 0, ie middle of screen)
    :param bounce_height: float - bounce height, see above 
        (defaults to value set in 'setup' routine's code component)
    :param period: float - seconds per complete bounce cycle
        (defaults to value set in 'setup' routine's code component)
    """
    # calculate new y coordinate
    # for info on `sin`, see comment in `pulsate_image` code above.
    new_y_coord = orig_y_coord + bounce_height * sin(time_val / period * 2 * pi)
    # assign new x-y position. note that this has to be passed as a 2-element
    # tuple, which is why we need to access the x-coordinate of the image component
    # by using `img.pos[0]` here.
    img.pos = (img.pos[0], new_y_coord)

# form a list of animation function names, which will be used
# to randomly decide, for each trial, which function/animation to use
att_animater_names = ['pulsate', 'flicker', 'bounce']
def get_dist_eyetracker(et):
    """
    A function for extracting participant('user') distance from 
    screen (not really the absolute distance, only the magnitude
    of the vector component that extends 'straight out of' the
    eyetracker, the distance in the 'y' axis according to Tobii). 
    Takes in an ioHub 'eyetracker device' object.
    Returns None if neither of the eyes were successfully captured.
    Returns a single numerical value representing mean distance of both eyes,
    or just distance of one if not both could be captured.
    """
    # `getLastSample` is not well documented, but includes basically all
    # data about last recording sample from eyetracker in a list.
    # the data that are relevant here happen to be at indices 15 and 34.
    # (at least for a Tobii Spectrum eyetracker - it's possible that
    # corresponding data will be found in other indices for other eyetrackers)
    last_sample = eye_tracker.getLastSample()
    # make extra sure that the last sample is valid, ie is a list of appropriate
    # length
    if not (last_sample and type(last_sample) == list and len(last_sample) > 34):
        return None
    # extract left/right eye Y axis distances from eyetracker
    left_dist = last_sample[15]
    right_dist = last_sample[34]
    left_truthy_not_nan = left_dist and not np.isnan(left_dist)
    right_truthy_not_nan = right_dist and not np.isnan(right_dist)
    if left_truthy_not_nan and right_truthy_not_nan:
        return (left_dist + right_dist) / 2
    if left_truthy_not_nan:
        return left_dist
    if right_truthy_not_nan:
        return right_dist
    return None



# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)

# Store info about the experiment session
psychopyVersion = '2021.2.3'
expName = 'audiovisual_volume_habituation_infant'  # from the Builder filename that created this script
expInfo = {'participant': '', 'session': '001'}
dlg = gui.DlgFromDict(dictionary=expInfo, sortKeys=False, title=expName)
if dlg.OK == False:
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
expInfo['psychopyVersion'] = psychopyVersion

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath='C:\\Users\\giobu365\\Documents\\psychopy_infant_audiovisual\\audiovisual_volume_habituation_infant.py',
    savePickle=True, saveWideText=True,
    dataFileName=filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.DEBUG)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp
frameTolerance = 0.001  # how close to onset before 'same' frame

# Start Code - component code to be run after the window creation

# Setup the Window
win = visual.Window(
    size=[1920, 1080], fullscr=True, screen=1, 
    winType='pyglet', allowGUI=False, allowStencil=False,
    monitor='spectrum_monitor', color=[1,1,1], colorSpace='rgb',
    blendMode='avg', useFBO=True, 
    units='deg')
# store frame rate of monitor if we can measure it
expInfo['frameRate'] = win.getActualFrameRate()
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess

# Setup eyetracking
ioDevice = ioConfig = ioSession = ioServer = eyetracker = None

# create a default keyboard (e.g. to check for escape)
defaultKeyboard = keyboard.Keyboard()

# Initialize components for Routine "setup"
setupClock = core.Clock()
# check if there is a pre-existing 'events.hdf5' file,
# meaning that a previous experiment run did
# not finish (since at the end of the experiment,
# the 'events.hdf5' file is renamed and moved to the
# data directory).
if os.path.isfile('events.hdf5'):
    # make sure that the 'incomplete_hdf5_files'
    # subdirectory in 'data' exists, and otherwise 
    # create it
    inc_dpath = os.path.join('data', 'incomplete_hdf5_files')
    if not os.path.isdir(inc_dpath):
        os.mkdir(inc_dpath)
    # generate a random identifier number/filename for
    # the file to be renamed with
    rand_num = randint(1, 10000000)
    inc_fname = f'inc_{rand_num}.hdf5'
    # rename the file and move it to
    # 'data/incomplete_hdf5_files' directory
    os.rename(
        'events.hdf5',
        os.path.join(inc_dpath, inc_fname)
     )

# specify name of YAML config file, located in the
# same directory as the .psyexp file.
# note that if 'mocking' is turned on ('Before
# Experiment' tab), the computer's mouse is
# used instead of an actual eyetracker.
if MOCK_ON:
    config_file = 'MOUSE_eyetracker.yaml'
else:
    config_file = 'tobii_config.yaml'

# specify name of eye tracker (ie the name
# specified in the YAML file)
eye_tracker_name = 'tracker'

# form file path to save eyetracker
# (ioHub/HDF5) data to
# (note that the file format, '.hdf5',
# is __not__ to be included in this file
# path, since ioHub automatically appends 
# '.hdf5' to the passed file path
new_hdf5_path = os.path.join(
        'data',
        '{}_{}_{}_hdf5'.format(expInfo['participant'],
            expInfo['expName'],
            expInfo['date'])
)

# attempt to connect to device specified in
# YAML config file
io_connection = launchHubServer(
    iohub_config_name=config_file,
    experiment_code='audiovisual_volume_habituation_infant',
    datastore_name=new_hdf5_path,
    window=win
)

# check if can get details for eye tracker device,
# and otherwise quit
if io_connection.getDevice(eye_tracker_name):
    # assign the tracker to a variable to make
    # it easier to reference
    eye_tracker = io_connection.getDevice(eye_tracker_name)
else:
    print(
        f"Could not connect to eye tracker '{eye_tracker_name}', is it on?"
        " Quitting..."
    )
    core.quit()


# Initialize components for Routine "start_screen"
start_screenClock = core.Clock()
text_start = visual.TextStim(win=win, name='text_start',
    text='Nu börjar det!',
    font='Open Sans',
    units='deg', pos=(0, 0), height=1.0, wrapWidth=None, ori=0.0, 
    color='black', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);

# Initialize components for Routine "trial"
trialClock = core.Clock()
text_trial1_timekeeper = visual.TextStim(win=win, name='text_trial1_timekeeper',
    text=None,
    font='Arial',
    pos=(0, 0), height=1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=0.0);
# draw a 'waiting screen' image+text and flip the window, to give the
# participant something to look at while stimuli are being loaded
load_screen_img = visual.ImageStim(
    win=win,
    name=f'load_screen_image', 
    image='stimuli/visual/load_screen_elephant.png', 
    mask=None,
    ori=0, 
    pos=(0, 0), 
    size=(LOAD_IMG_WIDTH, LOAD_IMG_HEIGHT),
    color=[1,1,1],
    colorSpace='rgb',
    opacity=1,
    flipHoriz=False, 
    flipVert=False,
    texRes=2048,
    interpolate=True,
    depth=0.0
)
load_screen_text = visual.TextStim(win=win, name='text_start',
    text=LOAD_MESSAGE,
    font='Open Sans',
    units='deg', pos=(0, -LOAD_IMG_HEIGHT/2 - MESSAGE_TXT_HEIGHT), 
    height=MESSAGE_TXT_HEIGHT, wrapWidth=None, ori=0.0, 
    color='black', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0
);

load_screen_img.draw()
load_screen_text.draw()
win.flip()

# -- pre-load social image stimuli by creating a list of
# image stimulus components --
# extract the subset of rows in the dataframe
# which pertain to social stimuli
social_img_df = img_df[img_df.type=='social']
# use function defined in the 'Before Experiment' tab to generate
# a list of PsychoPy image components based on the social images
img_ls_social = generate_img_components_from_df(
    social_img_df,
    img_category="social"
)

# pre-load non-social images, in the same way as was done for social images,
# repeating the steps for each type of non-social images
img_ls_geo = generate_img_components_from_df(
    img_df[img_df.type=='geometric'],
    img_category="geometric"
)
img_ls_manm = generate_img_components_from_df(
    img_df[img_df.type=='manmade'],
    img_category="manmade"
)
img_ls_nat = generate_img_components_from_df(
    img_df[img_df.type=='natural'],
    img_category="natural"
)


# all trial images should be at equal distances from
# center of screen, arranged in a 'cross-like'
# shape, ie they have equal x/y offsets from the
# vertical/horizontal center lines
x_offset = IMG_X_OFFSET
y_offset = IMG_Y_OFFSET
pos_ls = [
    (-x_offset, y_offset),
    (x_offset, y_offset),
    (-x_offset, -y_offset),
    (x_offset, -y_offset)
]

# a variable for a list of position indices, used for assigning the above 
# positions to trial visual stimuli. this list is updated for every trial,
# see the definition of 'generate_position_indices' in the 'Before
# Experiment' tab. the variable is initialized with a value of None 
# to ensure that an entirely random list of indices is generated for the
# first trial
position_indices = None

# pre-load attention grabber images by creating a list of
# attention grabber image components, specifying the
# attention image maximum width/height
att_img_ls = generate_img_components_from_df(
    att_img_df, 
    max_width=ATT_MAX_WIDTH, 
    max_height=ATT_MAX_WIDTH,
    with_lin_interpolation=False
)

# pre-load sound stimuli by creating a list of
# sound stimulus components
sound_ls = []
for row_ind in range(len(audio_df)):
    new_sound = sound.Sound(
        audio_df.iloc[row_ind]['file_path'],
        name=f'sound_{row_ind}', 
        stereo=True, 
        hamming=True,
        # set volume level to 0 as increasing from 0dB to 20dB then up to 70dB by 10dB
        # through 'assign_pseudorandom_volume' function
        volume=0
    )
    sound_ls.append(new_sound)

# shuffle, ie randomly reorder the list of sound components.
# it's used to go through all the sounds, playing each one across the different volume levels
# throught 7 trials in a row (including silent ones).
# Trial sequences for the different sounds are interleaved with the habituation task
shuffle(sound_ls)

# pre-load attention grabbing sounds by creating a list of
# sound stimulus components
att_sound_ls = []
for row_ind in range(len(att_audio_df)):
    new_att_sound = sound.Sound(
        att_audio_df.iloc[row_ind]['file_path'],
        name=f'att_sound_{row_ind}', 
        stereo=True,
        hamming=True,
        volume=1
    )
    att_sound_ls.append(new_att_sound)

# form list of trial 'area of interest components', 'invisible'
# rectangles that are used for assigning trial visual stimuli
# areas of interest, and checking if participant gaze
# is directed at/'located in' these AOI's
trial_aoi_components = []
for i in range(4):
    new_aoi_comp = visual.Rect(
        win=win,
        name=f'trial_aoi_comp_{i}', 
        ori=0, pos=(0, 0), 
        size=(TRIAL_AOI_WIDTH, TRIAL_AOI_HEIGHT),
        lineColor=[-1,-1,1], 
        colorSpace='rgb', 
        opacity=1,
        depth=0.0,
        interpolate=False
    )
    trial_aoi_components.append(new_aoi_comp)

# form attention grabber 'area of interest component'
att_aoi_comp = visual.Rect(
    win=win,
    name=f'att_aoi_comp', 
    ori=0, pos=(0, 0), 
    size=(ATT_AOI_WIDTH, ATT_AOI_HEIGHT),
    lineColor=[-1,1,-1], 
    colorSpace='rgb', 
    opacity=1,
    depth=0.0,
    interpolate=False
)
# since attention grabbers are always presented
# in the centre of the image, and have the same
# maximum width/height, the same AOI component
# can be used for all of them
for img in att_img_ls:
    img.aoi = att_aoi_comp

# there are issues in PsychoPy related to linear interpolation
# (ie when interpolate=True) for image components, which makes black lines
# appear around some of the images. 
# see eg
# https://discourse.psychopy.org/t/borders-on-some-images/9995/6
# for a discussion about this.
# the experiment however needs this interpolation to avoid having really 
# jagged image display.
# to get rid of the unwanted 'outlines' mentioned above, this white borders-only
# rectangle component is drawn on top of image components, to hide the outlines.
rect_outlinehider = visual.Rect(
    win=win,
    name='rect_outlinehider',
    units='deg',
    width=2,
    height=2,
    ori=0,
    pos=(0, 0),
    lineWidth=0.8,
    colorSpace='rgb',
    lineColor=[1, 1, 1],
    fillColor=None,
    opacity=None,
    depth=1,
    interpolate=False
)

# key response component for checking if experimenter hits
# 'a' key to play attention grabbing sound, or 'p' key to
# pause the experiment
key_resp_exp = keyboard.Keyboard()

# counter for number of trials having passed
trial_counter = 0

if DEBUG_ON:
    # text component for showing current trial time, for enabling
    # manual/external timing comparisons between eyetracker events
    # and events as seen on camera recording
    text_trialtime = visual.TextStim(win=win, name='text_trialtime',
        text='',
        font='Open Sans',
        units='deg', pos=(0, -6), height=0.8, wrapWidth=None, ori=0.0, 
        color='black', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=0.0
    )

# if eyetracker is being mocked (ie a mouse is used
# instead of an actual eyetracker), skip mirroring 
# display as described in comments below
if not MOCK_ON:
    # while running the experiment, the experimenter
    # is expected to keep an eye on where the participant
    # is looking and to play attention-grabbing sounds
    # (at beginning of trial, when attention grabber
    # image is displayed) if it seems necessary.
    # therefore, parts of the experiment need to be
    # mirrored to the experimenter's screen. hence,
    # a separate PsychoPy Window instance is created
    # here for the experiment's screen
    win_experimenter = visual.Window(
        size=[1920, 1080], 
        fullscr=True, 
        screen=0, 
        winType='pyglet', 
        allowGUI=False, 
        allowStencil=False,
        monitor='legion5_builtin', 
        color=[1,1,1], 
        colorSpace='rgb',
        blendMode='avg', 
        useFBO=True, 
        units='deg'
    )

    # a 'dot' for showing to the experimenter where
    # the participant is currently directing their
    # gaze
    gaze_dot = visual.GratingStim(
        win=win_experimenter,
        tex=None,
        mask='gauss',
        pos=(0, 0),
        size=(0.8, 0.8),
        color='green',
        units='deg'
    )

    # text stimuli for displaying participant
    # distance from tobii eyetracker
    part_dist_info_msg = visual.TextStim(
        win=win_experimenter,
        text="Participant distance from screen (target range 55-75):",
        pos=(0, -2),
        height=0.5,
        units='deg',
        autoLog=False,
        color='black'
    )
    part_dist_val_msg = visual.TextStim(
        win=win_experimenter,
        text="dist_val_msg",
        pos=(0, -3),
        height=0.5,
        units='deg',
        autoLog=False,
        color='black'
    )

    # message to show if participant distance from eyetracker could not
    # be recorded/measured by eyetracker
    missed_recording_str = 'Could not capture distance'


    # the following part section enables mirroring all trial stimuli to the experimenter screen
    # and is very technical (requires deeper understanding of Python). some more
    # info can be found at:
    # https://discourse.psychopy.org/t/a-general-approach-to-mirroring-all-stimuli-to-a-secondary-screen/22095

    # decorator function for modifying draw
    # methods so that they will draw to mirror (experimenter)
    # window as well
    def make_draw_mirror(draw_fun):
        def mirror_draw_fun(*args, **kwargs):
            draw_fun(win=win_experimenter)
            draw_fun(*args, **kwargs)
        return mirror_draw_fun

    # decorator function for making
    # a window flip method also trigger flipping of
    # mirror (experimenter) window
    def make_flip_mirror(flip_fun):
        def mirror_flip_fun(*args, **kwargs):
                win_experimenter.flip(*args, **kwargs)
                flip_fun(*args, **kwargs)
        return mirror_flip_fun

    # form list of all trial stimulus image components
    all_trial_img_ls = (
        att_img_ls +
        img_ls_social +
        img_ls_geo +
        img_ls_manm +
        img_ls_nat
    )

    for img_stim in all_trial_img_ls:
        # enforce drawing to both 'usual' and mirror
        # screen
        img_stim.draw = make_draw_mirror(img_stim.draw)


# Initialize components for Routine "trial_habit"
trial_habitClock = core.Clock()
text_trial2_timekeeper = visual.TextStim(win=win, name='text_trial2_timekeeper',
    text=None,
    font='Arial',
    pos=(0, 0), height=1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=0.0);
    
# Initialize components for Routine "end_screen"
end_screenClock = core.Clock()
text_end = visual.TextStim(win=win, name='text_end',
    text=END_MESSAGE,
    font='Open Sans',
    units='deg', pos=(0, 0), height=MESSAGE_TXT_HEIGHT, wrapWidth=None, ori=0.0, 
    color='black', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

# ------Prepare to start Routine "setup"-------
continueRoutine = True
# update component parameters for each repeat
# keep track of which components have finished
setupComponents = []
for thisComponent in setupComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
setupClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1

# -------Run Routine "setup"-------
while continueRoutine:
    # get current time
    t = setupClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=setupClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in setupComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "setup"-------
for thisComponent in setupComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# the Routine "setup" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# ------Prepare to start Routine "start_screen"-------
continueRoutine = True
routineTimer.add(2.000000)
# update component parameters for each repeat
# tell ioHub how to label the oncoming series of eye tracker
# data
io_connection.sendMessageEvent('start_recording')

# instruct eye tracker to stream eye data to PsychoPy's ioHub
# (meaning the data will be later stored in the participant's
# HDF5 data file)
eye_tracker.setRecordingState(True)

# keep track of which components have finished
start_screenComponents = [text_start]
for thisComponent in start_screenComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
start_screenClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1

# -------Run Routine "start_screen"-------
while continueRoutine and routineTimer.getTime() > 0:
    # get current time
    t = start_screenClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=start_screenClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *text_start* updates
    if text_start.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        text_start.frameNStart = frameN  # exact frame index
        text_start.tStart = t  # local t and not account for scr refresh
        text_start.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(text_start, 'tStartRefresh')  # time at next scr refresh
        text_start.setAutoDraw(True)
    if text_start.status == STARTED:
        # is it time to stop? (based on global clock, using actual start)
        if tThisFlipGlobal > text_start.tStartRefresh + 2.0-frameTolerance:
            # keep track of stop time/frame for later
            text_start.tStop = t  # not accounting for scr refresh
            text_start.frameNStop = frameN  # exact frame index
            win.timeOnFlip(text_start, 'tStopRefresh')  # time at next scr refresh
            text_start.setAutoDraw(False)
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in start_screenComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "start_screen"-------
for thisComponent in start_screenComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)

# set up handler to look after randomisation of conditions etc

trials_block_loop = data.TrialHandler(nReps=6.0, method='sequential', 
    extraInfo=expInfo, originPath=-1,
    trialList=[None],
    seed=None, name='trials_block_loop')
thisExp.addLoop(trials_block_loop)  # add the loop to the experiment
thisTrials_block_loop = trials_block_loop.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisTrials_pause_block_loop.rgb)
if thisTrials_block_loop != None:
    for paramName in thisTrials_block_loop:
        exec('{} = thisTrials_block_loop[paramName]'.format(paramName))

for thisTrials_block_loop in trials_block_loop:
    currentLoop = trials_block_loop
    # abbreviate parameter names if possible (e.g. rgb = thisTrials_block_loop.rgb)
    if thisTrials_block_loop != None:
        for paramName in thisTrials_block_loop:
            exec('{} = thisTrials_block_loop[paramName]'.format(paramName))
    
    # set up handler to look after randomisation of conditions etc
    trial_loop = data.TrialHandler(nReps=7.0, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=[None],
        seed=None, name='trial_loop')
    thisExp.addLoop(trial_loop)  # add the loop to the experiment
    thisTrial_loop = trial_loop.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisTrial_loop.rgb)
    if thisTrial_loop != None:
        for paramName in thisTrial_loop:
            exec('{} = thisTrial_loop[paramName]'.format(paramName))
    
    for thisTrial_loop in trial_loop:
        currentLoop = trial_loop
        # abbreviate parameter names if possible (e.g. rgb = thisTrial_loop.rgb)
        if thisTrial_loop != None:
            for paramName in thisTrial_loop:
                exec('{} = thisTrial_loop[paramName]'.format(paramName))
        
        # ------Prepare to start Routine "trial"-------
        continueRoutine = True
        # update component parameters for each repeat
        # send 'trial start' message to ioHub (ie the eye tracker), 
        # that will be included in the next eye tracker 
        # data record
        # NOTE that there are also some 'send message to ioHub' code
        # scattered throughout the 'code_trial' code snippet's
        # 'Each Frame' tab
        io_connection.sendMessageEvent(f'exp1 trial {trial_counter+1} start')
        
        # calculate start times/delays and durations of stimuli
        # attention grabber start time
        att_grab_extra_startdelay = random_20thsec((ATT_GRAB_LATEST - ATT_GRAB_EARLIEST)*20)
        att_grab_start = ATT_GRAB_EARLIEST + att_grab_extra_startdelay
        # auditory stimulus delay after gaze capture
        audio_extra_delay = random_20thsec((AUDITORY_MAX_DELAY - AUDITORY_MIN_DELAY)*20)
        audio_delay = AUDITORY_MIN_DELAY + audio_extra_delay
        # visual stimuli delay after audio onset, and visual stimuli duration
        visual_extra_delay = random_20thsec((VISUAL_MAX_DELAY - VISUAL_MIN_DELAY)*20)
        visual_delay = VISUAL_MIN_DELAY + visual_extra_delay
        visual_extra_dur = random_20thsec((VISUAL_MAX_DUR - VISUAL_MIN_DUR)*20)
        visual_dur = VISUAL_MIN_DUR + visual_extra_dur
        # duration of 'blank screen', shown after visual stimuli offset
        blank_extra_dur = random_20thsec((BLANK_MAX_DUR - BLANK_MIN_DUR)*20)
        blank_dur = BLANK_MIN_DUR + blank_extra_dur
        
        # booleans indicating if the different components are 
        # active/started (ie shown/played)
        att_active = False
        audio_started = False
        visual_started = False
        visual_active = False
        
        # list that represents the last 30 frames of the experiment, 
        # that's used for checking if participant's gaze has been directed
        # at attention grabber for 200ms out of the last 500ms.
        # 200ms corresponds to 12 frames, 500ms corresponds to 30 frames.
        # 0's in the list indicate that participant's gaze wasn't directed
        # at the attention grabber, while 1's mean that gaze _was_ directed
        # at attention grabber for a particular frame. 
        att_gaze_frame_ls = [0] * 30
        # boolean indicating if gaze has been directed at attention grabber 
        # for 200ms out of the last 500ms
        gaze_captured = False
        
        # timepoints (initialized with very high values,
        # to avoid events being preemptively triggered)
        # timepoint at which attention grabber starts
        # being presented
        att_onset_t = 9999999
        # timepoint at which gaze is directed at
        # attention grabber
        gaze_capture_t = 9999999
        # timepoint at which sound starts playing
        audio_onset_t = 9999999
        # timepoint at which visual stimuli start to be presented
        visual_onset_t = 9999999
        # timepoint at which visual stimuli stop being presented
        visual_end_t = 9999999
        
        # sounds are to be played at the same volume two trials in a row, meaning
        # trials come in pairs. accordingly, for all even-numbered trials 
        # (trial number 0, 2, 4...), a new sound stimulus,
        # and a new sound level, is to be used
        if trial_counter % 13 == 0:
            # determine which sound should be used, going from start to end of
            # list (it's 7 trials pop-out task for 7 volume levels + 6 trials habituation task = 13 trial complete block)
            sound_index = (trial_counter // 13) % len(sound_ls)
            sound_trial = sound_ls[sound_index]
            
        # if it's not a silent trial (in which case trial counter would be a multiple of 13 as first in the trial sequence) 
        # assign a volume pseudorandomly
        if 0 < trial_counter % 13 < 7:
            assign_volume(sound_trial)
        
        
        # generate a random array of visual stimuli to use
        img_ls = get_trial_img_ls()
        # assign positions to the visual stimuli, making sure that
        # each category's image is in a different position compared to where
        # the category's previous image (in the previous trial) appeared
        position_indices = generate_position_indices(position_indices)
        for image_index, position_index in enumerate(position_indices):
            img_ls[image_index].pos = pos_ls[position_index]
            # assign an 'area of interest' rectangle/component,
            # which is used for checking if gaze is directed at
            # the AOI around the stimulus' position
            img_ls[image_index].aoi = trial_aoi_components[position_index]
            img_ls[image_index].aoi.pos = pos_ls[position_index]
        
        # randomly pick an attention grabber image
        # to use for this trial
        att_img = randchoice(att_img_ls, 1)[0]
        # extract original size/opacity/y coordinate of attention grabber,
        # as they might change (depending on the type
        # of attention grabber animation used), meaning these original values
        # are necessary to reset the attention grabber component at 
        # end of routine/trial
        att_orig_size = att_img.size
        att_orig_opacity = att_img.opacity
        att_orig_y_coord = att_img.pos[1]
        
        # randomly pick which attention grabber
        # animation (function) to use
        animater_name = randchoice(att_animater_names, 1)[0]
        
        # reset variables/properties related to keyboard component 
        # (the one used for checking if
        # experimenter wants to play attention grabbing sound, or pause the
        # experiment)
        key_resp_exp.keys = []
        # boolean indicating that the experimenter pressed the 'a'
        # key last frame (set to True to prevent any previous routine's key press
        # from activating sound at first frame)
        soundpress_last_frame = True
        # variable for storing time at which last attention grabbing sound
        # started playing, to prevent experimenter from playing multiple 
        # attention grabbing sounds in a row. is initialized to a negative
        # value to ensure that no cooldown is applied at start.
        att_sound_start_time = -9999
        # variable for storing 'active' attention grabbing sound component,
        # for being able to stop it before trial start
        active_attention_sound = None
        # variable for counting the number of times that experimenter plays
        # an attention grabbing sound
        att_sound_counter = 0
        
        # dictionary for counting total number of frames for which visual stimuli
        # are shown, and how many of these frames are spent looking at each type
        # of visual stimulus
        gaze_counter_dict = {
            'total': 0,
            'manmade': 0,
            'natural': 0,
            'geometric': 0,
            'social': 0
        }
        
        # if debugging is enabled
        if DEBUG_ON:
            # boolean indicating if screenshot has been taken
            trial_screenshot_taken = False
            attgrab_screenshot_taken = False
        
        # get the time at which the first trial 'frame flip'
        # will occur, to save as trial start time
        trial_start_time = win.getFutureFlipTime(clock=None)
        
        # boolean indicating if experiment is currently paused (True) or not (False) 
        exp_paused = False
        # counter that indicates for how long this trial is paused during attention
        # grabbing phase. note that this will be a sum, if the trial is repeatedly
        # paused
        pause_duration = 0
        # boolean indicating if pause key was pressed down on last frame
        pausepress_last_frame = False
        # keep track of which components have finished
        trialComponents = [text_trial1_timekeeper]
        for thisComponent in trialComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        trialClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
        frameN = -1
        
        # -------Run Routine "trial"-------
        while continueRoutine:
            # get current time
            t = trialClock.getTime()
            tThisFlip = win.getFutureFlipTime(clock=trialClock)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *text_trial1_timekeeper* updates
            if text_trial1_timekeeper.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                text_trial1_timekeeper.frameNStart = frameN  # exact frame index
                text_trial1_timekeeper.tStart = t  # local t and not account for scr refresh
                text_trial1_timekeeper.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(text_trial1_timekeeper, 'tStartRefresh')  # time at next scr refresh
                text_trial1_timekeeper.setAutoDraw(True)
            # time to activate attention grabber?
            if t >= att_grab_start and not (att_active or gaze_captured):
                att_active = True
                att_t_start = tThisFlip
            
            # if experiment is paused, check if experimenter wants to unpause (otherwise
            # do nothing)
            if att_active and exp_paused:
                this_frame_keys = key_resp_exp.getKeys(
                    keyList=['a', 'A', 'p', 'P'],
                    waitRelease=False
                )
                # check if experimenter has pressed 'p' or 'P' to pause the experiment
                pausepress_this_frame = any([k.name.lower() == 'p' for k in this_frame_keys])
                new_pausepress = pausepress_this_frame and not pausepress_last_frame
                pausepress_last_frame = pausepress_this_frame
                # end pause
                if new_pausepress:
                    pause_duration += trialClock.getTime() - pause_trial_time
                    exp_paused = False
            
            # attention grabber active?
            if att_active and not exp_paused:
                # collect any experimenter keyboard input
                this_frame_keys = key_resp_exp.getKeys(
                    keyList=['a', 'A', 'p', 'P'],
                    waitRelease=False
                )
                # check if experimenter has pressed 'p' or 'P' to pause the experiment
                pausepress_this_frame = any([k.name.lower() == 'p' for k in this_frame_keys])
                new_pausepress = pausepress_this_frame and not pausepress_last_frame
                # update boolean so that it's possible to check if new
                # keyboard press occurs in next frame
                pausepress_last_frame = pausepress_this_frame
                if new_pausepress and not exp_paused:
                    # start pause
                    pause_trial_time = trialClock.getTime()
                    exp_paused = True
            #    if exp_paused:
            #        continue
                # time passed since attention grabber was activated
                att_passed_t = tThisFlip - att_t_start
                # animate the attention grabber image, based on which
                # animation function was randomly picked in 'Begin Routine'
                if animater_name == 'pulsate':
                    pulsate_image(
                        att_img,
                        att_passed_t,
                        att_orig_size
                    )
                elif animater_name == 'flicker':
                    flicker_image(
                        att_img,
                        att_passed_t,
                        att_orig_opacity
                    )
                elif animater_name == 'bounce':
                    bounce_image(
                        att_img,
                        att_passed_t,
                        att_orig_y_coord
                    )
                # show the attention grabber
                att_img.draw()
                # remove the 'oldest' frame's information about whether participant
                # gaze was then directed at attention grabber
                att_gaze_frame_ls.pop(0)
                # check if gaze (as measured by eye tracker)
                # is within attention grabber's area of interest (AOI)
                gaze_pos = eye_tracker.getPosition()
                if gaze_pos and att_img.aoi.contains(*gaze_pos, units='deg'):
                    # save information that participant gazed at attention grabber during
                    # this frame
                    att_gaze_frame_ls.append(1)
                else:
                    # save information that participant did _not_ gaze at
                    # attention grabber during this frame
                    att_gaze_frame_ls.append(0)
                # check for how many frames, out of the last 30, that participant gaze
                # was directed at attention grabber
                acc_gaze_frames = sum(att_gaze_frame_ls)
                # if participant gaze was directed at attention grabber for
                # at least 12 (200ms) out of the last 30 (500ms) frames, OR if maximum
                # attention grabber duration has passed, signal
                # that attention grabber phase is finished.
                if (acc_gaze_frames >= 12) or (att_passed_t - pause_duration) > ATT_GRAB_MAX_DUR:
                    att_active = False
                    gaze_captured = True
                    gaze_capture_t = tThisFlip
                    # send event marker to ioHub/'eyetracker output' that
                    # attention grabber phase is ending
                    io_connection.sendMessageEvent(f'exp1 trial {trial_counter+1} attention grabber end')
                # check if experimenter has pressed 'a' (also checking for caps 'A') 
                # key, meaning they want to play an attention grabbing sound
                soundpress_this_frame = any([k.name.lower() == 'a' for k in this_frame_keys])
                new_soundpress = soundpress_this_frame and not soundpress_last_frame
                # check if enough time has passed since last 
                # attention grabbing sound (if any)
                # started playing
                sound_cooldown_passed = (tThisFlip - att_sound_start_time) >= ATT_SOUND_COOLDOWN_TIME
                if new_soundpress and sound_cooldown_passed:
                    # make sure that if an attention grabbing sound has been played before,
                    # it is completely stopped before the new one is played
                    if active_attention_sound:
                        active_attention_sound.stop()
                    # randomly pick an attention grabbing sound to play
                    active_attention_sound = randchoice(att_sound_ls, 1, replace=False)[0]
                    active_attention_sound.play()
                    att_sound_start_time = tThisFlip
                    att_sound_counter += 1
                # update boolean so that it's possible to check if new
                # keyboard press occurs in next frame
                soundpress_last_frame = soundpress_this_frame
            
            # gaze has been captured and time of delay
            # from gaze capture->sound onset has passed?
            if gaze_captured and not audio_started and t >= (gaze_capture_t + audio_delay):
                # if an attention grabbing sound is playing, stop it
                if active_attention_sound:
                    active_attention_sound.stop()
                audio_started = True
                # if it's not a silent trial, play the trial's sound component/audio
                # stimulus
                if sound_trial is not None:
                    sound_trial.play()
                audio_onset_t = tThisFlip
                # send event marker to ioHub/'eyetracker output' that
                # trial sound is starting (or would have started, if it weren't
                # a silent trial)
                io_connection.sendMessageEvent(f'exp1 trial {trial_counter+1} sound onset')
            
            # gaze has been captured and time of delay
            # from sound onset->visual stimul onset has passed?
            if audio_started and not visual_started and tThisFlip >= (audio_onset_t + visual_delay):
                visual_active = True
                visual_started = True
                visual_onset_t = tThisFlip
                # send event marker to ioHub/'eyetracker output' that
                # trial visual array is shown (on the next frame flip)
                io_connection.sendMessageEvent(f'exp1 trial {trial_counter+1} visual onset')
            
            # visual stimuli are active?
            if visual_active:
                # increment counter for number of frames that visual stimuli are
                # presented for
                gaze_counter_dict['total'] += 1
                for img in img_ls:
                    # draw the image itself
                    img.draw()
                    # copy the image's width/height to the 'outline hider'
                    # rectangle component, and draw it to hide
                    # unwanted and interpolation-related outlines
                    # (this is a hack to circumvent a bug - see the 'begin routine'
                    # tab)
                    rect_outlinehider.width = img.size[0]
                    rect_outlinehider.height = img.size[1]
                    rect_outlinehider.pos = img.pos
                    rect_outlinehider.draw()
                # check if gaze (as measured by eye tracker)
                # is within any visual stimulus' area of interest (AOI)
                gaze_pos = eye_tracker.getPosition()
                if gaze_pos:
                    for img in img_ls:
                        if img.aoi.contains(*gaze_pos, units='deg'):
                            # increase counter for number of frames
                            # spent by participant gazing at image
                            # category's images
                            gaze_counter_dict[img.img_category] += 1
                # visual stimuli have been shown for
                # set duration?
                if (tThisFlip - visual_onset_t) >= visual_dur:
                    visual_active = False
                    visual_end_t = tThisFlip
                    # send event marker to ioHub/'eyetracker output' that
                    # trial visual array is ending (on the next frame flip)
                    io_connection.sendMessageEvent(f'exp1 trial {trial_counter+1} visual offset')
            
            # end blank screen has been shown for its
            # intended duration?
            if tThisFlip >= (visual_end_t + blank_dur):
                # end trial
                continueRoutine = False
            
            # if debugging is enabled
            if DEBUG_ON:
                # draw areas of interest outlines
                att_aoi_comp.draw()
                for img in img_ls:
                    img.aoi.draw()
                # NOTE that taking screenshots of trials, as below, causes very noticeable
                # hiccups in experiment flow/screen refresh time. if this is a problem
                # for timing measurements, comment out the screenshot-related code
                # take screenshot of trial at attention grabber onset
                if not attgrab_screenshot_taken and att_active:
                    win.getMovieFrame(buffer='back')
                    win.saveMovieFrames(f'exp_snapshots/att_grab_{trial_counter}.png', codec='png')
                    attgrab_screenshot_taken = True
                # take screenshot of trial at visual stimuli onset
                if not trial_screenshot_taken and visual_active:
                    win.getMovieFrame(buffer='back')
                    win.saveMovieFrames(f'exp_snapshots/trial_{trial_counter}.png', codec='png')
                    trial_screenshot_taken = True
                # show current experiment/PsychoPy time
                text_trialtime.text = tThisFlipGlobal
                text_trialtime.draw()
            
            # if eyetracker is being mocked (ie a mouse is used
            # instead of an actual eyetracker), skip mirroring 
            # display as described in comments
            if not MOCK_ON:
                # draw areas of interest outlines
                att_aoi_comp.draw(win=win_experimenter)
                for img in img_ls:
                    img.aoi.draw(win=win_experimenter)
            
                # get latest gaze position coordinates
                gaze_pos_formirror = eye_tracker.getPosition()
            
                # check if a valid gaze position could be
                # recorded (in which case the returned value is
                # of type tuple or list, as this is how the
                # coordinates are 'packaged')
                gaze_pos_isvalid = isinstance(
                    gaze_pos_formirror,
                    (tuple, list)
                )
                
                # update and show message
                # which indicates participant distance from
                # eyetracker
                y_dist = get_dist_eyetracker(eye_tracker)
                if y_dist:
                    part_dist_val_msg.setText(y_dist)
                else:
                    # if a value is shown currently, replace it with
                    # indication that distance could not be measured
                    if part_dist_val_msg.text != missed_recording_str:
                        part_dist_val_msg.setText(missed_recording_str)
                part_dist_info_msg.draw()
                part_dist_val_msg.draw()
            
                # update marker dot to match gaze position
                if gaze_pos_isvalid:
                    gaze_dot.pos = gaze_pos_formirror
                gaze_dot.draw()
            
                # now that everything has been drawn on the
                # window's 'canvas' , 'flip' the window/screen
                # so that the new information is displayed
                # to the experimenter
                win_experimenter.flip()
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in trialComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # -------Ending Routine "trial"-------
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # send 'trial end' message to ioHub/the eye tracker 
        # that will be included in the next eye tracker 
        # data record
        io_connection.sendMessageEvent(f'exp1 trial {trial_counter+1} end')
        
        # get trial end time
        trial_end_t = t
        
        # assign attention grabber image component its 
        # original size/opacity/position, in case they 
        # were manipulated for animation (depends on what
        # type of animation was used)
        att_img.size = att_orig_size
        att_img.opacity = att_orig_opacity
        att_img.pos = (att_img.pos[0], att_orig_y_coord)
        
        ### BEGIN SAVE TRIAL DATA ###
        # trial start time (counting from start of experiment)
        thisExp.addData('trial_global_start_time', trial_start_time)
        
        # the time information below is relative to experiment start
        
        # intended (ie as instructed to PsychoPy - should be in agreement
        # with recorded actual times) start times/durations/delays
        # (intended) attention grabber start time
        thisExp.addData('att_grab_start_time_intended', att_grab_start)
        # (intended) auditory stimulus delay after gaze capture
        thisExp.addData('gaze_to_audio_delay_intended', audio_delay)
        # (intended) visual stimuli delay after audio onset
        thisExp.addData('audio_to_visual_delay_intended', visual_delay)
        # (intended) visual stimuli duration
        thisExp.addData('visual_duration_intended', visual_dur)
        # (intended) duration of 'blank screen', shown after visual stimuli offset
        thisExp.addData('end_blank_duration_intended', blank_dur)
        
        # actual (recorded during experiment) time information
        # (actual) attention grabber start time
        thisExp.addData('att_grab_start_time_actual', att_t_start)
        # (actual) gaze captured time (ie time at which participant has
        # looked at attention grabber for 200 out of the last 500ms, and
        # attention grabber is removed)
        thisExp.addData('gaze_captured_time', gaze_capture_t)
        # (actual) auditory stimulus onset time
        thisExp.addData('audio_onset_time', audio_onset_t)
        # (actual) visual stimuli onset time
        thisExp.addData('visual_onset_time', visual_onset_t)
        # (actual) visual stimuli offset time
        thisExp.addData('visual_offset_time', visual_end_t)
        # (actual) trial end time
        thisExp.addData('trial_end_time', trial_end_t)
        
        # number of times that experimenter played an attention grabbing
        # sound
        thisExp.addData('attention_sounds_played', att_sound_counter)
        
        # total number of frames during which stimuli were shown
        # (there are normally 60 frames per second)
        thisExp.addData('visual_stimuli_duration_nframes', gaze_counter_dict['total'])
        # proportion of visual stimuli presentation time spent by participant
        # gazing at the stimuli AOI's, rounded to five decimals
        thisExp.addData(
            'visual_social_prop', 
            round(gaze_counter_dict['social'] / gaze_counter_dict['total'], 5)
        )
        thisExp.addData(
            'visual_geometric_prop', 
            round(gaze_counter_dict['geometric'] / gaze_counter_dict['total'], 5)
        )
        thisExp.addData(
            'visual_manmade_prop', 
            round(gaze_counter_dict['manmade'] / gaze_counter_dict['total'], 5)
        )
        thisExp.addData(
            'visual_natural_prop', 
            round(gaze_counter_dict['natural'] / gaze_counter_dict['total'], 5)
        )
        
        # visual stimuli characteristics
        for img in img_ls:
            # visual stimuli file paths
            thisExp.addData(
                f'visual_{img.img_category}_filepath', img.image
            )
            # visual stimuli position x/y coordinates, in degrees
            thisExp.addData(
                f'visual_{img.img_category}_pos_x', img.pos[0]
            )
            thisExp.addData(
                f'visual_{img.img_category}_pos_y', img.pos[1]
            )
        
        # audio stimuli characteristics
        # if it was a silent trial, save 'silent' as sound filepath,
        # and '0' for volume
        if sound_trial is None:
            thisExp.addData(
                f'audio_filepath', 'silent',
            )
            thisExp.addData(
                f'audio_volume', 0,
            )
        else:
            thisExp.addData(
                f'audio_filepath', sound_trial.sound,
            )
            thisExp.addData(
                f'audio_volume', sound_trial.volume,
            )
        
        # total (attention grabbing phase) pause duration
        thisExp.addData('pause_duration', pause_duration)
        
        ### END SAVE TRIAL DATA ###
        
        # increment trial counter
        trial_counter += 1
        
        # the Routine "trial" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        thisExp.nextEntry()
        
    # completed 7.0 repeats of 'trial_loop' along the different volume levels
    # get the randomized video path
    video_df = pd.read_csv(HVIDEO_SPEC_PATH)
    # 'in-trial' audio stimuli
    haudio_df = pd.read_csv(HAUDIO_SPEC_PATH)

    if habit_trial_list[a] % 2 == 0: 
            HVIDEO_FILE = video_df[video_df.type == 'crash']
            fileindx = 0

    if habit_trial_list[a] % 2 == 1: 
        HVIDEO_FILE = video_df[video_df.type == 'nocrash']
        fileindx = 1

    
    # set up handler to look after randomisation of conditions etc
    trial_habit_loop = data.TrialHandler(nReps=6.0, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=[None],
        seed=None, name='trial_habit_loop')
    thisExp.addLoop(trial_habit_loop)  # add the loop to the experiment
    thisTrial_habit_loop = trial_habit_loop.trialList[0]  # so we can initialise stimuli with some values
    
    # abbreviate parameter names if possible (e.g. rgb = thisTrial_loop.rgb)
    if thisTrial_habit_loop != None:
        for paramName in thisTrial_habit_loop:
            exec('{} = thisTrial_habit_loop[paramName]'.format(paramName))
        
    for thisTrial_habit_loop in trial_habit_loop:
        currentLoop = trial_habit_loop
        # abbreviate parameter names if possible (e.g. rgb = thisTrial_loop.rgb)
        if thisTrial_habit_loop != None:
            for paramName in thisTrial_habit_loop:
                exec('{} = thisTrial_habit_loop[paramName]'.format(paramName))
        
        # ------Prepare to start Routine "trial_habit"-------
        continueRoutine = True
        routineTimer.add(9.000000) # removed as did not continue
        
        # update component parameters for each repeat
        movie_habituation = visual.MovieStim3(
            win=win, name='movie_habituation',units='deg', 
            noAudio = True,
            filename=HVIDEO_FILE.file_path[fileindx],
            ori=0.0, pos=(0, 0), opacity=None,
            loop=False,
            size=(HVIDEO_WIDTH_DEG, HVIDEO_HEIGHT_DEG),
            depth=0.0,
            )
        sound_habit_trial = sound.Sound(
            haudio_df.iloc[fileindx]['file_path'],
            name=f'sound_{fileindx}', 
            stereo=True, 
            hamming=True,
            volume=1
        )    
        # send 'trial start' message to ioHub (ie the eye tracker), 
        # that will be included in the next eye tracker 
        # data record
        # NOTE that there are also some 'send message to ioHub' code
        # scattered throughout the 'code_trial' code snippept's
        # 'Each Frame' tab
        io_connection.sendMessageEvent(f'exp2 habituation trial {trial_counter+1} start')
        
        # calculate start times/delays and durations of stimuli
        # attention grabber start time
        att_grab_extra_startdelay = random_20thsec((ATT_GRAB_LATEST - ATT_GRAB_EARLIEST)*20)
        att_grab_start = ATT_GRAB_EARLIEST + att_grab_extra_startdelay
        
        # auditory stimulus delay after gaze capture
        audio_extra_delay_habit = random_20thsec((AUDITORY_MAX_DELAY_HABIT - AUDITORY_MIN_DELAY_HABIT)*20)
        audio_delay_habit = AUDITORY_MIN_DELAY_HABIT + audio_extra_delay_habit
        
        # visual stimuli delay after gaze capture
        video_extra_delay = random_20thsec((VIDEO_MAX_DELAY - VIDEO_MIN_DELAY)*20)
        video_delay = VIDEO_MIN_DELAY + video_extra_delay
        
        video_extra_dur = random_20thsec((VIDEO_MAX_DUR - VIDEO_MIN_DUR)*20)
        video_dur = VIDEO_MIN_DUR + video_extra_dur
        
        # duration of 'blank screen', shown after visual stimuli offset
        blank_extra_dur = random_20thsec((BLANK_MAX_DUR - BLANK_MIN_DUR)*20)
        blank_dur = BLANK_MIN_DUR + blank_extra_dur
        
        # booleans indicating if the different components are 
        # active/started (ie shown/played)
        att_active = False
        audio_started = False
        video_started = False
        video_active = False
        
        # list that represents the last 30 frames of the experiment, 
        # that's used for checking if participant's gaze has been directed
        # at attention grabber for 200ms out of the last 500ms.
        # 200ms corresponds to 12 frames, 500ms corresponds to 30 frames.
        # 0's in the list indicate that participant's gaze wasn't directed
        # at the attention grabber, while 1's mean that gaze _was_ directed
        # at attention grabber for a particular frame. 
        att_gaze_frame_ls = [0] * 30
        # boolean indicating if gaze has been directed at attention grabber 
        # for 200ms out of the last 500ms
        gaze_captured = False
        
        # timepoints (initialized with very high values,
        # to avoid events being preemptively triggered)
        # timepoint at which attention grabber starts
        # being presented
        att_onset_t = 9999999
        # timepoint at which gaze is directed at
        # attention grabber
        gaze_capture_t = 9999999
         # timepoint at which sound starts playing
        audio_onset_t = 9999999
        # timepoint at which visual stimuli start to be presented
        video_onset_t = 9999999
        # timepoint at which visual stimuli stop being presented
        video_end_t = 9999999
        
        # randomly pick an attention grabber image
        # to use for this trial
        att_img = randchoice(att_img_ls, 1)[0]
        # extract original size/opacity/y coordinate of attention grabber,
        # as they might change (depending on the type
        # of attention grabber animation used), meaning these original values
        # are necessary to reset the attention grabber component at 
        # end of routine/trial
        att_orig_size = att_img.size
        att_orig_opacity = att_img.opacity
        att_orig_y_coord = att_img.pos[1]
        
        # randomly pick which attention grabber
        # animation (function) to use
        animater_name = randchoice(att_animater_names, 1)[0]
        
        # reset variables/properties related to keyboard component 
        # (the one used for checking if
        # experimenter wants to play attention grabbing sound, or pause the
        # experiment)
        key_resp_exp.keys = []
        # boolean indicating that the experimenter pressed the 'a'
        # key last frame (set to True to prevent any previous routine's key press
        # from activating sound at first frame)
        soundpress_last_frame = True
        # variable for storing time at which last attention grabbing sound
        # started playing, to prevent experimenter from playing multiple 
        # attention grabbing sounds in a row. is initialized to a negative
        # value to ensure that no cooldown is applied at start.
        att_sound_start_time = -9999
        # variable for storing 'active' attention grabbing sound component,
        # for being able to stop it before trial start
        active_attention_sound = None
        # variable for counting the number of times that experimenter plays
        # an attention grabbing sound
        att_sound_counter = 0
        
       
        # if debugging is enabled
        if DEBUG_ON:
            # boolean indicating if screenshot has been taken
            trial_screenshot_taken = False
            attgrab_screenshot_taken = False
        
        # get the time at which the first trial 'frame flip'
        # will occur, to save as trial start time
        trial_habit_start_time = win.getFutureFlipTime(clock=None)
        
        # boolean indicating if experiment is currently paused (True) or not (False) 
        exp_paused = False
        # counter that indicates for how long this trial is paused during attention
        # grabbing phase. note that this will be a sum, if the trial is repeatedly
        # paused
        pause_duration = 0
        # boolean indicating if pause key was pressed down on last frame
        pausepress_last_frame = False
        
        # keep track of which components have finished
        habituation_trialComponents = [text_trial2_timekeeper]
        for thisComponent in habituation_trialComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
                
        habituation_videoComponents = [movie_habituation]
        for thisComponent in habituation_videoComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        trial_habitClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
        frameN = -1

        # -------Run Routine "trial_habit"-------
        while continueRoutine and routineTimer.getTime() > 0:
            # get current time
            tThisFlip = win.getFutureFlipTime(clock=trial_habitClock)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
                        
            # *text_trial2_timekeeper* updates
            if text_trial2_timekeeper.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                text_trial2_timekeeper.frameNStart = frameN  # exact frame index
                text_trial2_timekeeper.tStart = t  # local t and not account for scr refresh
                text_trial2_timekeeper.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(text_trial2_timekeeper, 'tStartRefresh')  # time at next scr refresh
                text_trial2_timekeeper.setAutoDraw(True)
                
            # time to activate attention grabber?
            if t >= att_grab_start and not (att_active or gaze_captured):
                att_active = True
                att_t_start = tThisFlip
            
            # if experiment is paused, check if experimenter wants to unpause (otherwise
            # do nothing)
            if att_active and exp_paused:
                this_frame_keys = key_resp_exp.getKeys(
                    keyList=['a', 'A', 'p', 'P'],
                    waitRelease=False
                )
                # check if experimenter has pressed 'p' or 'P' to pause the experiment
                pausepress_this_frame = any([k.name.lower() == 'p' for k in this_frame_keys])
                new_pausepress = pausepress_this_frame and not pausepress_last_frame
                pausepress_last_frame = pausepress_this_frame
                # end pause
                if new_pausepress:
                    pause_duration += trialClock.getTime() - pause_trial_time
                    exp_paused = False
            
            # attention grabber active?
            if att_active and not exp_paused:
                # collect any experimenter keyboard input
                this_frame_keys = key_resp_exp.getKeys(
                    keyList=['a', 'A', 'p', 'P'],
                    waitRelease=False
                )
                # check if experimenter has pressed 'p' or 'P' to pause the experiment
                pausepress_this_frame = any([k.name.lower() == 'p' for k in this_frame_keys])
                new_pausepress = pausepress_this_frame and not pausepress_last_frame
                # update boolean so that it's possible to check if new
                # keyboard press occurs in next frame
                pausepress_last_frame = pausepress_this_frame
                if new_pausepress and not exp_paused:
                    # start pause
                    pause_trial_time = trialClock.getTime()
                    exp_paused = True
            #    if exp_paused:
            #        continue
                # time passed since attention grabber was activated
                att_passed_t = tThisFlip - att_t_start
                # animate the attention grabber image, based on which
                # animation function was randomly picked in 'Begin Routine'
                if animater_name == 'pulsate':
                    pulsate_image(
                        att_img,
                        att_passed_t,
                        att_orig_size
                    )
                elif animater_name == 'flicker':
                    flicker_image(
                        att_img,
                        att_passed_t,
                        att_orig_opacity
                    )
                elif animater_name == 'bounce':
                    bounce_image(
                        att_img,
                        att_passed_t,
                        att_orig_y_coord
                    )
                # show the attention grabber
                att_img.draw()
                
                #win.flip() ### CAN BE REMOVED WHEN NOT MOCKED (???)
                
                # remove the 'oldest' frame's information about whether participant
                # gaze was then directed at attention grabber
                att_gaze_frame_ls.pop(0)
                # check if gaze (as measured by eye tracker)
                # is within attention grabber's area of interest (AOI)
                gaze_pos = eye_tracker.getPosition()
                if gaze_pos and att_img.aoi.contains(*gaze_pos, units='deg'):
                    # save information that participant gazed at attention grabber during
                    # this frame
                    att_gaze_frame_ls.append(1)
                else:
                    # save information that participant did _not_ gaze at
                    # attention grabber during this frame
                    att_gaze_frame_ls.append(0)
                # check for how many frames, out of the last 30, that participant gaze
                # was directed at attention grabber
                acc_gaze_frames = sum(att_gaze_frame_ls)
                # if participant gaze was directed at attention grabber for
                # at least 12 (200ms) out of the last 30 (500ms) frames, OR if maximum
                # attention grabber duration has passed, signal
                # that attention grabber phase is finished.
                if (acc_gaze_frames >= 12) or (att_passed_t - pause_duration) > ATT_GRAB_MAX_DUR:
                    att_active = False
                    gaze_captured = True
                    gaze_capture_t = tThisFlip
                    # send event marker to ioHub/'eyetracker output' that
                    # attention grabber phase is ending
                    io_connection.sendMessageEvent(f'exp2 habituation trial {trial_counter+1} attention grabber end')
                # check if experimenter has pressed 'a' (also checking for caps 'A') 
                # key, meaning they want to play an attention grabbing sound
                soundpress_this_frame = any([k.name.lower() == 'a' for k in this_frame_keys])
                new_soundpress = soundpress_this_frame and not soundpress_last_frame
                # check if enough time has passed since last 
                # attention grabbing sound (if any)
                # started playing
                sound_cooldown_passed = (tThisFlip - att_sound_start_time) >= ATT_SOUND_COOLDOWN_TIME
                if new_soundpress and sound_cooldown_passed:
                    # make sure that if an attention grabbing sound has been played before,
                    # it is completely stopped before the new one is played
                    if active_attention_sound:
                        active_attention_sound.stop()
                    # randomly pick an attention grabbing sound to play
                    active_attention_sound = randchoice(att_sound_ls, 1, replace=False)[0]
                    active_attention_sound.play()
                    att_sound_start_time = tThisFlip
                    att_sound_counter += 1
                # update boolean so that it's possible to check if new
                # keyboard press occurs in next frame
                soundpress_last_frame = soundpress_this_frame
            
            
            # gaze has been captured and time of delay
            # from gaze capture->sound onset has passed?
            if gaze_captured and movie_habituation.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                
                # if an attention grabbing sound is playing, stop it
                if active_attention_sound:
                    active_attention_sound.stop()
                video_started = True
                video_active = True
                video_onset_t = tThisFlip
                
                # keep track of start time/frame for later
                movie_habituation.frameNStart = frameN  # exact frame index
                movie_habituation.tStart = t  # local t and not account for scr refresh
                movie_habituation.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(movie_habituation, 'tStartRefresh')  # time at next scr refresh
                movie_habituation.setAutoDraw(True)
                
                # send event marker to ioHub/'eyetracker output' that
                # trial sound is starting (or would have started, if it weren't
                # a silent trial)
                io_connection.sendMessageEvent(f'exp2 trial {trial_counter+1} video onset')
            
                           
            # visual stimuli are active?
            if video_active:
                
                # gaze has been captured and time of delay
                # from gaze capture->sound onset has passed?
                if not audio_started and (tThisFlip - video_onset_t) >= audio_delay_habit:
                    audio_started = True
                    sound_habit_trial.play()
                    audio_onset_t = tThisFlip
                    # send event marker to ioHub/'eyetracker output' that
                    # trial sound is starting (or would have started, if it weren't
                    # a silent trial)
                    io_connection.sendMessageEvent(f'exp2 trial {trial_counter+1} sound onset')
                
                # visual stimuli have been shown for
                # set duration?
                if (tThisFlip - video_onset_t) >= video_dur:
                    video_active = False
                    video_end_t = tThisFlip
                    movie_habituation.tStop = t  # not accounting for scr refresh
                    movie_habituation.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(movie_habituation, 'tStopRefresh')  # time at next scr refresh
                    movie_habituation.setAutoDraw(False)
                    # send event marker to ioHub/'eyetracker output' that
                    # trial visual array is ending (on the next frame flip)
                    io_connection.sendMessageEvent(f'exp2 trial {trial_counter+1} video offset')
            if movie_habituation.status == FINISHED:  # force-end the routine
                continueRoutine = False
                        
            # end blank screen has been shown for its
            # intended duration?
            if tThisFlip >= (video_end_t + blank_dur):
                # end trial
                continueRoutine = False
            
            # if eyetracker is being mocked (ie a mouse is used
            # instead of an actual eyetracker), skip mirroring 
            # display as described in comments
            if not MOCK_ON:
                
                # get latest gaze position coordinates
                gaze_pos_formirror = eye_tracker.getPosition()
            
                # check if a valid gaze position could be
                # recorded (in which case the returned value is
                # of type tuple or list, as this is how the
                # coordinates are 'packaged')
                gaze_pos_isvalid = isinstance(
                    gaze_pos_formirror,
                    (tuple, list)
                )
                
                # update and show message
                # which indicates participant distance from
                # eyetracker
                y_dist = get_dist_eyetracker(eye_tracker)
                if y_dist:
                    part_dist_val_msg.setText(y_dist)
                else:
                    # if a value is shown currently, replace it with
                    # indication that distance could not be measured
                    if part_dist_val_msg.text != missed_recording_str:
                        part_dist_val_msg.setText(missed_recording_str)
                part_dist_info_msg.draw()
                part_dist_val_msg.draw()
            
                # update marker dot to match gaze position
                if gaze_pos_isvalid:
                    gaze_dot.pos = gaze_pos_formirror
                gaze_dot.draw()
            
                # now that everything has been drawn on the
                # window's 'canvas' , 'flip' the window/screen
                # so that the new information is displayed
                # to the experimenter
                win_experimenter.flip()
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                break
            continueRoutine = False  # will revert to True if at least one component still running
            
            for thisComponent in habituation_trialComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            for thisComponent in habituation_videoComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # -------Ending Routine "trial_habit"-------
        for thisComponent in habituation_trialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
                
        for thisComponent in habituation_videoComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        movie_habituation.stop()
       
        # send 'habituation trial end' message to ioHub/the eye tracker 
        # that will be included in the next eye tracker 
        # data record
        io_connection.sendMessageEvent(f'exp2 habituation trial {trial_counter+1} end')
        
        # get trial end time
        habit_trial_end_t = t
        
        # assign attention grabber image component its 
        # original size/opacity/position, in case they 
        # were manipulated for animation (depends on what
        # type of animation was used)
        att_img.size = att_orig_size
        att_img.opacity = att_orig_opacity
        att_img.pos = (att_img.pos[0], att_orig_y_coord)
        
        ### BEGIN SAVE TRIAL DATA ###
        # trial start time (counting from start of experiment)
        thisExp.addData('trial_habituation_global_start_time', trial_habit_start_time)
        
        # the time information below is relative to experiment start
        
        # intended (ie as instructed to PsychoPy - should be in agreement
        # with recorded actual times) start times/durations/delays
        # (intended) attention grabber start time
        thisExp.addData('att_grab_start_time_intended', att_grab_start)
        # (intended) auditory stimulus delay after gaze capture
        thisExp.addData('gaze_to_audio_delay_intended', audio_delay_habit)
        # (intended) visual stimuli delay after gaze capture
        thisExp.addData('video_to_attention_delay_intended', video_delay)
        # (intended) duration of 'blank screen', shown after visual stimuli offset
        thisExp.addData('end_blank_duration_intended', blank_dur)
        
        # actual (recorded during experiment) time information
        # (actual) attention grabber start time
        thisExp.addData('att_grab_start_time_actual', att_t_start)
        # (actual) gaze captured time (ie time at which participant has
        # looked at attention grabber for 200 out of the last 500ms, and
        # attention grabber is removed)
        thisExp.addData('gaze_captured_time', gaze_capture_t)
        # (actual) video stimuli onset time
        thisExp.addData('video_onset_time', video_onset_t)
        # (actual) auditory stimulus onset time
        thisExp.addData('audio_onset_time', audio_onset_t)
        # (actual) visual stimuli offset time
        thisExp.addData('video_offset_time', video_end_t)
        # (actual) trial end time
        thisExp.addData('trial_end_time', habit_trial_end_t)
        
        # number of times that experimenter played an attention grabbing
        # sound
        thisExp.addData('attention_sounds_played_habituation', att_sound_counter)
        
        # habituation video
        thisExp.addData('habituation_video_fpath', HVIDEO_FILE.file_path[fileindx])
        thisExp.addData(f'habituation_video_type', HVIDEO_FILE.type[fileindx])
        
        # audio stimuli characteristics
        # if it was a silent trial, save 'silent' as sound filepath,
        # and '0' for volume
        thisExp.addData(f'habituation_audio_filepath', sound_habit_trial.sound)
        thisExp.addData(f'habituation_audio_volume', sound_habit_trial.volume)
            
        # total (attention grabbing phase) pause duration
        thisExp.addData('pause_duration', pause_duration)
        
        ### END SAVE TRIAL DATA ###
        
        # increment trial counter
        trial_counter += 1
        
        # the Routine "trial" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        thisExp.nextEntry()
    a += 1      
    # completed 6.0 repeats of 'trial_habit_loop' over habituation repetitions

    # completed 6.0 repeats of 'trial_block_loop' over block repetitions

# ------Prepare to start Routine "end_screen"-------
continueRoutine = True
routineTimer.add(5.000000)
# update component parameters for each repeat
# keep track of which components have finished
end_screenComponents = [text_end]
for thisComponent in end_screenComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
end_screenClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1

# -------Run Routine "end_screen"-------
while continueRoutine and routineTimer.getTime() > 0:
    # get current time
    t = end_screenClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=end_screenClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *text_end* updates
    if text_end.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        text_end.frameNStart = frameN  # exact frame index
        text_end.tStart = t  # local t and not account for scr refresh
        text_end.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(text_end, 'tStartRefresh')  # time at next scr refresh
        text_end.setAutoDraw(True)
    if text_end.status == STARTED:
        # is it time to stop? (based on global clock, using actual start)
        if tThisFlipGlobal > text_end.tStartRefresh + 5.0-frameTolerance:
            # keep track of stop time/frame for later
            text_end.tStop = t  # not accounting for scr refresh
            text_end.frameNStop = frameN  # exact frame index
            win.timeOnFlip(text_end, 'tStopRefresh')  # time at next scr refresh
            text_end.setAutoDraw(False)
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in end_screenComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "end_screen"-------
for thisComponent in end_screenComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
thisExp.addData('text_end.started', text_end.tStartRefresh)
thisExp.addData('text_end.stopped', text_end.tStopRefresh)
# stop eye tracker recording
eye_tracker.setRecordingState(False)
# close io connection
io_connection.quit()

# Flip one final time so any remaining win.callOnFlip() 
# and win.timeOnFlip() tasks get executed before quitting
win.flip()

# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename+'.csv', delim='auto')
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
