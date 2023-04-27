---
layout: page
title: Week 1A Intro to Shell - Navigating Files and Directories
---

Navigating Files and Directories
================================

> Overview
> --------
>
> 
> **Questions**
> 
> *   How can I perform operations on files outside of my working directory?
>     
> *   What are some navigational shortcuts I can use to make my work more efficient?
>     
> 
> **Objectives**
> 
> *   Use a single command to navigate multiple steps in your directory structure, including moving backwards (one level up).
>     
> *   Perform operations on files in directories outside your working directory.
>     
> *   Work with hidden directories and hidden files.
>     
> *   Interconvert between absolute and relative paths.
>     
> *   Employ navigational shortcuts to move around your file system.
>     

Moving around the file system
-----------------------------

We’ve learned how to use `pwd` to find our current location within our file system. We’ve also learned how to use `cd` to change locations and `ls` to list the contents of a directory. Now we’re going to learn some additional commands for moving around within our file system.

Use the commands we’ve learned so far to navigate to the `data` directory, if you’re not already there.

    $ cd [your scratch directory]
    $ cd data
    

What if we want to move back up and out of this directory and to our top level directory? Can we type `cd shell_data`? Try it and see what happens.

    $ cd shell_data
    
    -bash: cd: shell_data: No such file or directory
    

Your computer looked for a directory or file called `shell_data` within the directory you were already in. It didn’t know you wanted to look at a directory level above the one you were located in.

We have a special command to tell the computer to move us back or up one directory level.

    $ cd ..
    

Now we can use `pwd` to make sure that we are in the directory we intended to navigate to, and `ls` to check that the contents of the directory are correct.

    $ pwd    

    $ ls
     

From this output, we can see that `..` did indeed take us back one level in our file system.

You can chain these together like so:

    $ ls ../../
    
- What does this show?

![Permissions breakdown](../assets/img/filesystem-challenge.svg)


> Finding hidden directories
> --------------------------
> 
> First navigate to your home directory. There is a hidden directory within this directory. Explore the options for `ls` to find out how to see hidden directories. List the contents of the directory and identify the name of the text file in that directory.
> 
> Hint: hidden files and folders in Unix start with `.`, for example `.my_hidden_directory`
> 
> > Solution
> > --------
> > 
> > First use the `man` command to look at the options for `ls`.
> > 
> >     $ man ls
> >     
> > 
> > The `-a` option is short for `all` and says that it causes `ls` to “not ignore entries starting with .” This is the option we want.
> > 
> >     $ ls -a
> >     
> >     
> > 
> > And then list the contents of the directory using `ls`.
> > 
> >     $ ls
> >     
> >     
> > 

In most commands the flags can be combined together in no particular order to obtain the desired results/output.

    $ ls -Fa
    $ ls -laF
    

Examining the contents of other directories
-------------------------------------------

By default, the `ls` commands lists the contents of the working directory (i.e. the directory you are in). You can always find the directory you are in using the `pwd` command. However, you can also give `ls` the names of other directories to view. Navigate to your home directory if you are not already there.

    $ cd
    

Full vs. Relative Paths
-----------------------

The `cd` command takes an argument which is a directory name. Directories can be specified using either a _relative_ path or a full _absolute_ path. The directories on the computer are arranged into a hierarchy. The full path tells you where a directory is in that hierarchy. Navigate to the home directory, then enter the `pwd` command.

    $ cd data
    $ pwd  
    

Now enter the following command:

    $ cd [scratch location]/data
    

These two commands have the same effect, they both take us to the `data` directory. The first uses the absolute path, giving the full address from the home directory. The second uses a relative path, giving only the address from the working directory. A full path always starts with a `/`. A relative path does not.

A relative path is like getting directions from someone on the street. They tell you to “go right at the stop sign, and then turn left on Main Street”. That works great if you’re standing there together, but not so well if you’re trying to tell someone how to get there from another country. A full path is like GPS coordinates. It tells you exactly where something is no matter where you are right now.

You can usually use either a full path or a relative path depending on what is most convenient. If we are in the home directory, it is more convenient to enter the full path. If we are in the working directory, it is more convenient to enter the relative path since it involves less typing.

Over time, it will become easier for you to keep a mental note of the structure of the directories that you are using and how to quickly navigate amongst them.



> Key Points
> ----------
> 
> *   The `/`, `.`, and `..` characters represent important navigational shortcuts.
>     
> *   Hidden files and directories start with `.` and can be viewed using `ls -a`.
>     
> *   Relative paths specify a location starting from the current location, while absolute paths specify a location from the root of the file system.
>     

* * *

Adapted from the Data Carpentry Intro to Command Line -shell genomics https://datacarpentry.org/shell-genomics/
