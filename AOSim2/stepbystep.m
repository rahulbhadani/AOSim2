% A simple script set up as a step by step guide to AOSim2
% 
% Author: Alexander Rodack
% Date: 2/16/2015
% Software provided by Johanan L. Codona
clear all; clc; close all;
%% Getting Started!
% It is important to keep the updated version of AOSim2 as the version you
% are using.  The easiest way to do this is to make use of command line git
% commands.  GUIs also exist if that floats your boat, but they can be
% somewhat of a crutch.  It is always good to learn a new skill, and using
% command line is a powerful way to do things.
%
% This is easier on Linux than on Windows, and I have no Apple computers,
% but I assume it would be closer to Linux than Windows.  This will provide
% some simple instructions based on Linux, but all the git commands will
% work no matter where you use them as long as you have git installed on
% your computer.
%
% Open the command line tool (Terminal in Linux).
% If you don't have git, type: sudo apt-get install git
% If you have git, you are ready to make a clone of the AOSim2/OPTI528
% branch.
% To do this, simply type: git clone -b OPTI528 https://github.com/jlcodona/AOSim2.git
% This will create a directory with the name AOSim2. If you want to specify
% the name of the directory and where it is located, add a path to the end
% of the command above.  An example is: /home/alex/Desktop/AOSim_git/
% After you type in that command, press enter, and you will some lines
% appearing in Terminal.  There shouldn't be any errors, and the last line
% will read "Checking connectivity...done."
%
% BOOM! You know have AOSim2 on your computer
% 
% To examine this a bit, use the cd command to switch into the directory
% that was just created.  Now type: git branch
% You should see "* OPTI528"  This means everything worked correctly.
% Now type: git status
% You will see: On branch OPTI528. Your branch is up-to-date with 'origin/OPTI528'.
% This means you have the most recent version of the repo.
%
% Now that you have the software, we will examine how to use it.  There
% will be more on git later.


%% Object Oriented Matlab
% For those of you who are not familiar with Object Oriented Programming
% (OOP), this will provide a little insight.  Mathworks has some good
% tutorial videos (I know from experience) that will walk you through the
% basics of creating a class structure, and using it. I highly recommend
% them if you have never seen OOP before.
%
% Basically, in OOP, you create objects that are of user defined class
% structures.  These are essentially data structures that know how to use
% the data within them to do something.  
%
% This comes to fruition in the classes by means of properties and methods.
% The properties are the data items that store whatever you want (scalars,
% vectors, matrices, cells, strings, doubles, logicals, EVERYTHING that
% is/was/can be done in Matlab).  The methods are the functions that tell
% the class what to do with that stored data.
%
% If you look at a class in AOSim2, or in the tutorial videos, you will see
% that it is broken up into sections.  It will start with the name of the
% class and the superclass structure it is going to inherit its initial
% properties and methods from.  This can be handle, matlab.mixin.Copyable,
% or the name of any user defined class (you will see this a lot in AOSim2)
% Next is a list of properties. There are protected and public properties.
% If you want to learn more, look into the tutorials.  Following the
% properties list are the methods.  The first method is always called the
% Constructor. This method tells the class how to create an object of its
% class type.  This is the place to start when looking at a new class so
% you can understand how to create an object to use in your scripting.  All
% the methods that follow (sometimes called utilities) add functionality to
% the class.  This is the part that makes the data structure know what to
% do with the data stored within it.  Finally, there are overloaded methods
% (methods that have the same name as a method in a superclass to change
% them to fit the current class) and static methods. Again, if you want to
% know more, do the tutorials. They are helpful. Really.
%
% One more note, which you would know if you watched the tutorials or are
% already in the know about OOP. You call a property or method by using a
% dot operator.  This is why you will see a lot of g.name or g.coords.
% Basically, if the property or method is known to the class, you can
% access it by using the dot.

%% Let's Jump In! Creating the (AO)Grid
% Alright, now that some of the logistical stuff is out of the way, lets do
% some stuff.
%
% AOGrid is located in the @AOGrid folder in the AOSim2 section of the
% repository.  The reaon for all the '@' symbols in front of the folder
% names is that it allows for other functions to be included in the folder
% and act as methods in the class, but also maintain the ability to be
% called separate from the class.
%
% Now that you have AOGrid open, take a moment to read the comments at the
% beginning.  These are some notes from John that are important to using
% all of AOSim2 correctly.  Also, John has some great comments, so read
% them.  Some are helpful in understanding what is happening, others are
% good for a laugh or two.
%
% Now look over the properties.  This are things that will be a part of
% every class in AOSim2 (every other class has AOGrid as a superclass).
% Some are fairly straightforward, others might be somewhat harder to
% understand what they are for.  Knowing the properties of AOGrid is less
% important than for some of the other classes, but you should get a basic
% understanding of them.  We will learn more about them as we run into the
% need to use them.  Now take a look at the Constructor method.  You can
% see here how to create an AOGrid object, and the possible input arguments
% that are accepted by the code.  Finally, take some time to look at all
% the other methods that are included.  There again is no real reason to
% explain them all here, and a lot of them are used when others are called.
% I will go into detail about some of them when the need arises.
%
% So let's make a simple AOGrid:
Grid = AOGrid(64);

% This creates an AOGrid object that contains a 256x256 array.  The default
% values for data properties seen in the Constructor method are set to the
% object as well. 

% Let's explore a couple of the methods that can quickly become important
% coords and COORDS
% Lets run them and see what we get:
[x,y] = Grid.coords;
[X,Y] = Grid.COORDS;

% Take a look at the results
% you will see that x,y are vectors, and X,Y are meshgrid-like matrices.
% These commands map a real life coordinate system (in meters) to the
% pixels in the AOGrid based on the spacing property (defaults to 0.04)
% and how many pixels you gave to the array when creating the object (64
% in this example script).
% This are very convenient for things like plotting, and for understanding
% what a simulation might be like in physical units.

% grid
% This is another important method. This is used a lot in other classes
% that come down the line, but is useful to look at here.  Calling
% Grid.grid will print the current array stored in the object. If you have
% a matrix already created somewhere that is the right size, you can set
% the property grid_ (the actual array is stored here) to that matrix.  You
% can also set grid_ by using the constant method.

Grid.grid
matrix_A = magic(64);
Grid.grid(matrix_A);
Grid.grid
Grid.constant(1);
Grid.grid