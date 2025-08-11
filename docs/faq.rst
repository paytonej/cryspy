What is the _a pyori_ project philosophy?
-----------------------------------------

We intend for _a pyori_ to be a cross-platform set of tools for EBSD data analysis that is implemented primarily in Python packages that are installed with Python(x,y) and Spyder.

All functions in the OV library are directly available for the user. There are very few private functions. You can mix and match the functions that you need in order to do any analysis that you can dream up.


Why Python?
-----------

It doesn't take very long to pick up the basics of playing guitar and be able to play a tune, yet even seasoned musicians continue to use guitars to make music. Python is similar: it is a relatively easy language to pick up for people with little programming experience, and professional programmers routinely turn to Python to make sophisticated programs. With packages like python(x,y), inexperienced users have an extremely powerful environment at their fingertips to try out new ideas. The environment is free and extensible, and some knowledge of Python will remain useful to the beginner.

The primary purpose of the orientation visualizer library is to allow for new EBSD data analysis techniques to be invented and implemented by researchers interested in doing such things. A secondary purpose is to allow for automated analysis in routine investigations through scripting.


What is the difference between _a pyori_ and mtex?
-------------------------------------------

The core focus of MTEX is the mathematics of texture analysis, while the core focus of _a pyori_ is enabling the testing of novel EBSD data analysis ideas. MTEX requires a MATLAB license, while _a pyori_ is written in Python, which is both free in both the 'beer' and 'speech' senses.

Mtex offers many functions that are not currently available in OV (but may be available in the future). It is important to note that OV does not currently support texture analysis, 3D data sets, non-cubic pole figures, tensor analyses, and many other things that are available in mtex.

What is the difference between _a pyori_ and DREAM.3D?
------------------------------------------------------

The core focus of DREAM.3D is the analysis of 3D EBSD datasets, while _a pyori_ only has the capability of working with 2D data.

Why can't the best features of all of these EBSD freeware projects be combined?
-------------------------------------------------------------------------------

OV is essentially a library of functions that are available to the user, with a simple user interface that allows for routine analyses. 

When developing new analysis techniques or implementing experimental ones, it is desirable to have full access to a library of subfunctions. The OV user interface is built off of this library. The OV library is written with the purpose of making each function as intuitive as possible, and most of the code is extensively commented so that interested users can easily understand the code.

OV intentionally avoids using user-defined object classes, so the user has the possibility of manipulating their data as they wish. This is incompatible with mtex, since object classes are a fundamental part of mtex. It would take a lot of work to make every aspect of every function and class in mtex fully available for the user – this would require a complete rewriting of mtex. On the other hand, implementing mtex algorithms into OV is relatively easy.

One major drawback of using specially-defined object classes is that the end user needs to have a greater understanding of the operations of the program.

Mtex is covered by the GPLv2. This is a copyleft license that allows the code to be reused as long as the end program is also free and open source. Therefore, code can be shared between OV and mtex.

Disclaimers
-----------

OV has been developed in the course of the regular work duties of its authors. This code is a work in progress and is distributed on an as-is basis without warranties or conditions of any kind, either express or implied.