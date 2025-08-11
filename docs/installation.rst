These instructions explain how to get OVLib up and running on your computer. In principle, any distribution of Python 2.7 will do; however, the instructions below are written specifically for using the Spyder IDE and Python(x,y). We recommend the Spyder IDE due to its Matlab-like interface, which many users may find familiar.

1. Install Python(x,y) and visvis
=================================

Detailed instructions for Windows users
---------------------------------------
1. Download and install `python(x,y) <https://code.google.com/p/pythonxy/>`_
2. Install the `visvis plugin <http://code.google.com/p/pythonxy/wiki/AdditionalPlugins>`_ for Python(x,y).
3. Install `pyqtgraph <http://www.pyqtgraph.org/>` using the Windows installer.

Detailed instructions for Mac users
-----------------------------------
Installing python(x,y) on a macintosh requires compiling it from the source code. This is *very* easy; however, downloading and compiling all the components can take several hours to complete.

Follow these steps:

1. Install Apple's Xcode Developer Tools
	a. Launch Xcode and accept the EULA
	b. If OS is older than Mountain Lion, you may need to go into preferences and add command line tools.
2. Install MacPorts using the .pkg installer.
3. Open up a terminal window and run the following
	a. sudo port selfupdate
	b. sudo port upgrade outdated
	c. sudo port -v install py27-spyder
	d. sudo port install py27-opengl
	e. sudo easy-install pip
	f. pip install visvis
	g. pip install pyqtgraph
4. Drag the spyder.app (available from me) to your applications folder, double click to run it.

2. Install OVLib using Mercurial
================================

1. Install a Mercurial client

* For Mac users, we recommend `MacHg <https://bitbucket.org/jfh/machg/wiki/Home>`_
* For Windows users, we recommend: `TortoiseHg <http://tortoisehg.bitbucket.org/>`_

2. Clone the OVLib bitbucket repository. The following links may be helpful:

* http://ratfactor.com/mercurial-named-branches#bitbucket
* http://sce.uhcl.edu/support/files/StepsMacHG.pdf

To update your OVLib, pull incoming changes.

Contributing to ovlib may be done by creating a bitbucket account, forking the code, making your own changes, then pushing them to bitbucket.

3. Add OVLib to your Python Path
================================
1. Open up Spyder, go to the "Tools" menu, and click on PYTHONPATH manager.
2. Add your local ovlib directory to your python path and click "sychronize."
3. Restart Spyder


Beginner's tutorial
===================
*Coming soon.*