Development Process
==========================

Git and GitHub
^^^^^^^^^^^^^^^^^
bioScience is hosted on GitHub, and in order to contribute, you must create a `free GitHub account <https://github.com/signup/free>`_. `Git <https://git-scm.com/>`_ is used for version control to enable multiple people to collaborate on a project. If you are unfamiliar with Git, you can consult the `Git documentation <https://git-scm.com/doc>`_ to learn more.

In addition, the project adheres to a forking protocol that is described in greater detail on this page, in which contributors fork the repository, make modifications, and then create a pull request. Please read and adhere to all the instructions in this manual.

If you are unfamiliar with contributing to projects on GitHub through forking, consult the GitHub documentation for contributing to projects. GitHub provides a short tutorial using a test repository that may help you become more familiar with forking a repository, cloning a fork, creating a feature branch, pushing changes and making pull requests. The following resources can assist you in learning more about forking and pull requests on GitHub:

* the `GitHub documentation for forking a repo <https://docs.github.com/en/get-started/quickstart/fork-a-repo>`_.
* the `documentation on GitHub for collaborating on merge requests <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests>`_.
* the `documentation for working with variants on GitHub <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks>`_.

There are `instructions on GitHub <https://docs.github.com/en/get-started/quickstart/set-up-git>`_ for installing git, configuring SSH keys, and configuring git. Before you can work seamlessly between your local repository and GitHub, you must perform these procedures.

Create a fork of bioScience
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To modify the code, you will need your own copy of bioScience (also known as a fork). Visit the page for the `bioScience project <https://github.com/aureliolfdez/bioscience>`_ and click the ``Fork`` button. Before selecting ``Create Fork``, please deselect the option to reproduce only the main branch. You will clone your machine's fork.

.. code-block:: bash

    git clone https://github.com/your-user-name/bioscience.git bioscience-yourname
    cd bioscience-yourname
    git remote add upstream https://github.com/aureliolfdez/bioscience.git
    git fetch upstream

This creates the directory ``bioscience-yourname`` and connects your repository to the upstream (main project) bioScience repository.

Create a feature branch on GitHub
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Your local ``main`` branch should always reflect the current state of bioScience repository. First ensure it's up-to-date with the main bioScience repository.

.. code-block:: bash

    git checkout main
    git pull upstream main --ff-only

Then, create a feature branch for making your changes. For example:

.. code-block:: bash

    git checkout -b my-feature-branch

This changes your working branch from ``main```to the ``my-feature-branch`` branch. Keep any changes in this branch specific to one bug or new feature so it is clear what the branch brings to library. You can have many feature branches and switch in between them using the ``git checkout`` command.

When you want to update the feature branch with changes in main after you created the branch, check the section :any:`updatingPull`.

Making code changes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For files you intended to modify or add, run.

.. code-block:: bash

    git add <file>

Once you have made code changes, you can see all the changes you’ve currently made by running.

.. code-block:: bash

    git status

Finally, commit your changes to your local repository with an explanatory commit message

.. code-block:: bash

    git commit -m "commit message"

Pushing your code changes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
When you want your changes to appear publicly on your GitHub page, push your forked feature branch’s commits:

.. code-block:: bash

    git push origin my-feature-branch

Here ``origin`` is the default name given to your remote repository on GitHub. You can see the remote repositories

.. code-block:: bash

    git remote -v

If you added the upstream repository as described above you will see something like

.. code-block:: bash

    origin  git@github.com:yourname/bioscience.git (fetch)
    origin  git@github.com:yourname/bioscience.git (push)
    upstream git://github.com/aureliolfdez/bioscience.git (fetch)
    upstream git://github.com/aureliolfdez/bioscience.git (push)


.. _makingPull:

Making a pull request on GitHub
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
One you have finished your code changes, you are ready to make a pull request. A pull request is how code from your local repository becomes available to the GitHub community to review and merged into project to appear the in the next release. To submit a pull request:

1. Navigate to your repository on GitHub.
2. Click on the ``Compare & pull request`` button.
3. You can then click on ``Commits`` and ``Files Changed`` to make sure everything looks okay one last time.
4. Write a descriptive title that includes prefixes. bioScience uses a convention for title prefixes. Here are some common ones along with general guidelines for when to use them.
5. Write a description of your changes in the ``Preview Discussion`` tab.
6. Click on the ``Create pull request`` button.

This request then goes to the repository maintainers, and they will review the code.

.. note::
    Prefixes to write a descriptive title:
    
    * **BUG:** Bug fix.
    * **BLD:** Updates to the build process/scripts.
    * **DOC:** Additions/updates to documentation.
    * **FUN:** New functionality.
    * **HPC:** Performance improvement (High-Performance Computing)
    * **PER:** Performance issues. 
    * **QUE:** Submit question. 


.. _updatingPull:

Updating your pull request
^^^^^^^^^^^^^^^^^^^^^^^^^^
It is also important that updates in the bioScience main branch are reflected in your pull request. To update your feature branch with changes in the bioScience main branch, run:

.. code-block:: bash

    git checkout my-feature-branch
    git fetch upstream
    git merge upstream/main

If there are no conflicts (or they could be fixed automatically), a file with a default commit message will open, and you can simply save and quit this file.

If there are merge conflicts, you need to solve those conflicts. See for example at https://help.github.com/articles/resolving-a-merge-conflict-using-the-command-line/ for an explanation on how to do this.

Once the conflicts are resolved, run:

1. ``git add -u`` to stage any files you've updated;
2. ``git commit`` to finish the merge.

After the feature branch has been update locally, you can now update your pull request by pushing to the branch on GitHub:

.. code-block:: bash

    git push origin my-feature-branch

Any ``git push`` will automatically update your pull request with your branch's changes.

Updating the development environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
It is essential to regularly update your local main branch with updates from the bioScience main branch and to update your development environment to reflect any changes to the packages used during development.

.. code-block:: bash
    
    git checkout main
    git merge upstream/main
    # activate the virtual environment based on your platform
    python -m pip install --upgrade -r requirements.txt


.. _successfullPull:

Tips for a successfull pull request
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you have made it to the :any:`makingPull`, the core team may take a look. To improve the chances of your pull request being reviewed, you should:

* Ensure you have appropriate tests. These should be the first part of any pull request.
* Keep your pull requests as simple as possible. Larger pull requests take longer to review.
* Ensure that continuous integration is in a green state. Reviewers may not even look otherwise.
* Keep :any:`updatingPull`, either by request or every few days.

