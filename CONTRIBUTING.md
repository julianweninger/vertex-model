# Utopia Contribution and Support Guidelines

This document is both a guideline for how you can contribute to the Utopia project, but also of what is expected of Utopia users to get full support for Utopia.
It aims to streamline the process of using and developing in and for Utopia and keep us connected as a group, working together on this project.

To that end, there are a number of tasks that should be carried out periodically, e.g. every two weeks.
Every line marked with a ❇️ denotes such a task.


### Keep your Utopia installation and branches up to date
- ❇️ [<= 2 weeks] Pull the latest master in the [`utopia/utopia`](https://ts-gitlab.iup.uni-heidelberg.de/utopia/utopia) repository

    ```bash
    cd utopia
    git checkout master && git pull
    cd build && cmake ..
    ```

- ❇️ [<= 2 weeks] Pull the latest master in the [`utopia/vertex-model`](https://ts-gitlab.iup.uni-heidelberg.de/utopia/vertex-model) repository

    ```bash
    cd vertex-model
    git checkout master && git pull
    cd build && cmake ..
    ```
    
- ❇️ [<= 2 weeks] Merge the latest master into your model branch

    ```bash
    git checkout my-model-branch
    git merge master
    git push
    ```
    
- ❇️ [~ daily] Push relevant changes

    ```bash
    git checkout my-model-branch
    git push
    ```

_Note:_ Staying up-to-date with the latest master is usually not a big deal and we advise to do this as regularly as possible; in fact, the further you diverge from the master, the harder the update becomes.  
The [Utopia Mattermost channel](https://ts-chat.iup.uni-heidelberg.de/all/channels/utopia) informs about new features in the `utopia/utopia` master.
For the `utopia/vertex-model` repo, you can see in your model's merge request how many commits you are "behind master".


### Have a version of your model in the master
❇️ [~ monthly] Merge a stable version of your model into the `utopia/vertex-model` master.

During the first weeks of development, you might not deem your model stable.
However, there usually comes a point where you don't make any _major_ changes any more; that is a good point to pause, do some clean-up, write tests, get feedback, and eventually merge the model into the master.

**Benefits** once your model is in the master:

- It is tested alongside framework changes; thereby, it is made sure that changes to the framework code do not break your model
- You get feedback on your implementation
- You have a much cleaner and more focussed version control work flow as your following changes can be more granular

**To get your model merged:**

1. Open a Merge Request (if you don't already have one)
1. Fill in the MR description, making use of the existing templates
    * Write a one-sentence summary of what this MR does.
    * Fill in the checklist; you may delete points that you don't need or ignore them if you are unclear about what they mean.
1. The pipeline needs to be passing. This means:
    * Your model builds without warnings
    * Your model tests (if you implemented them) pass
    * If both these points are fulfilled, there will be only green ticks in the MR description, yay :tada:
1. Once ready, remove the `WIP:` in the title, and assign someone for review
    * You can assign basically anyone to review it
    * If you want feedback on something specific, you can comment directly on the changed code
    * If you want your model merged _without_ an in-depth review, denote that in the MR and assign @herdeanu, @lriedel, @hmack or @yunus.
      They will then merely assert that the implementation fulfils basic requirements, e.g. does not break any other parts of the repository.


### Contribute to Utopia or other vertex-model
- ❇️ [Every two weeks] Review a merge request
    - Depending on your level of expertise and your interests, this can be a model review or some framework MR.
    - Don't be shy; everyone starts somewhere and every input can be valuable.
    - Benefits for you: learn about clever implementations, best practices, testing methodologies, ... 
    - If you like, do it more frequently :)
- ❇️ [Whenever observed] [Report](https://ts-gitlab.iup.uni-heidelberg.de/utopia/utopia/issues/new?issue) errors, bugs, or ideas for improvements
- ❇️ [Whenever a feature is needed] [Write an issue](https://ts-gitlab.iup.uni-heidelberg.de/utopia/utopia/issues/new?issue)
    - If confident, propose an implementation by opening a merge request
    - You will get assistance... and lots of thanks from other users that also desire such a feature! :)

