**Welcome to the hack week!**

Login to datalore to access the jupyter notebook -- http://datalore.io

Login to slack to access the workspace --  http://slack.com

The slack has the overleaf link under the `#paper` channel.

--------------------------
Just in case it's helpful, a (very) quick guide to the basics of using github:

**Getting this repository on your computer**  
I) Create a github account, log in  
II) Return to the gradhackweek repo  
Hit the green button above: 'Clone or download' --> 'Clone with HTTPS' --> copy link  
III) Copy the repo to your computer,

```
cd <parent directory for your local copy of the repo>   
git clone <copied link>  
cd gradhackweek/
```

**Creating and working on a new branch**  
You're now on the master (main) branch. Branches are copies of the repo that can  
exist with different content in parallel (e.g., so we can work on different parts of the  
project simultaneously).

If you don't yet have a separate branch you're working on, create one,

`git checkout -b <branch name>`

Here 'checkout' moves you to <branch name>. '-b' creates a new branch.

Now keep your files for the project in this directory.

**Pushing local changes to the github repo**  
Push them back to the github repo as often as needed (at least at the end of each day during the hack  
week so we can integrate changes). You can push back to github in a few steps:

I) Show the state of your local repo compared to the branch on github,

`git status`

II) If you're happy with the local changes, add them for staging to the commit you're about to make,

`git add .`

Here '.' is to add all files in local directory; alternatively replace '.' with <filename>

III) Record the changes in a commit (keep the quotes around the commit message),

`git commit -m <"Informative commit message">`

IV) Push the commit to the github repo,

`git push origin <branch name>`

**Downloading and optionally merging the github repo with your local copy**  
To update your local repo with the current version on github, the safest way is to first just  
download the remote version of your branch without merging them into your local branch,

`git fetch`

This downloads all the remote branches; alternatively git fetch origin <branch name>

Compare your local repo to the verison you just downloaded

`git diff`

It's safest to save a version of your local copy before merging the downloaded version with it,

`git stash`

Now if you want, merge the downloaded version of the branch with your local one,

`git merge`

That's it!
