''' kevent.py

 Used to modify matplotlib plots to allow quitting simply by pressing
 the 'w' key for a single window or 'q' to close all windows
'''

import matplotlib.pyplot as window

# -----------------------------------------------------------------------------
def press(event):
    ''' Quit if q is pressed, close window if w is pressed. '''
    if event.key == 'q':
		window.close('all')
    
    if event.key == 'w': 
		window.close()
