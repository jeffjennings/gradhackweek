Welcome to the hack week!

The slack is here -- http://bit.ly/gradhackslack
    The slack has logistical info (participant list, daily schedule, etc.)
    
The overleaf is here -- http://bit.ly/gradhackms

--------------------------
Just in case it's helpful, a (very) quick guide to the basics of using github:

1) Create a github account, log in
2) Return to the gradhackweek repository.
Hit the green button above: 'Clone or download' --> 'Clone with HTTPS' --> copy link
3) Copy the repo to your computer:

$ cd <parent directory for your local copy of the repo>
$ git clone <copied link>
$ cd gradhackweek/

Awesome, you're now on the master (main) branch. Branches are copies of the repo that can 
exist with different content in parallel (e.g., so we can work on different parts of the 
project simultaneously). 

If you don't yet have a separate branch you're working on, create one,

$ git checkout -b <branch name> # 'checkout' moves you to <branch name>. '-b' creates a new branch

Now keep your files for the project in this directory, and push them back to the github repo 
as often as needed (at least at the end of each day during the hack week so we can 
integrate changes). You can push back to github in a few steps:

1) Show the state of your local repo compared to the branch on github,

$ git status 

2) If you're happy with the local changes, add them for staging to the commit you're about to make,

$ git add . # '.' to add all files in local directory; alternatively replace '.' with <filename>

3) Record the changes in a commit (keep the quotes around the commit message), 

$ git commit -m "Informative commit message" 

4) Push the commit to the github repo,

$ git push origin <branch name>

That's it!
