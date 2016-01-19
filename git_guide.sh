###	TUTORIAL	###
https://www.atlassian.com/git/tutorials/setting-up-a-repository
https://git-scm.com/book/en/v2/Getting-Started-About-Version-Control





#	SET-UP	#

# transform the current directory into a git repository
git init

# create a central repository where changes can't be made
git init --bare

# clone an existing repository
git clone <repo> # tends to be a github link,  defaults remote branch name as 'origin/master'

# define the authorname, author email & alias for log
git config --global user.name <name>
git config --global user.email <email>
git config --global alias.graph 'log --graph --oneline --full-history --all --decorate'
# git graph
git config --global alias.list '!ls'		# git list = ls, use ! for external commands



#	COMMITS	#

# select changes to files to be recorded/saved within git
git add <file/dir>	# select/stage; track and staged to be commited
git commit 		# record # git message format: summary of changes in <50 chars then dashed lines to detail changes
git commit --amend	# combine new changes with the most recent commit
git commit *		# everything in dir except hiddenfiles

# create a list of files to ignore when adding or commiting
touch .gitignore # use text editor i.e. *pyc

# see changes that have taken place since last commit
git status
git status -s	# short version

# remove file from working dir and from git tracking
rm <file>; git rm <file>	# requires a git commit
git rm --staged <file>		# remove from staging area

# change file name
mv <file> <file_renamed>
git mv <file> <file_renamed>

# display commited snapshots
git log --stat				# view changes between commits (-p flag also)
git log -p -2				# -p like -stat but more detail, -2 last 2 commits only
git log --oneline --decorate		# quickly summarise commits
git log --oneline master..<branch2>	# compare two branches, showing where difference is in commits
git log --graph --oneline --full-history --all --decorated # tree showing all branches
git log --pretty=format:"%h - %an, %ar : %s"	# show log with certain info only
git log --since=1.day
git log --grep delete.em		# filter only commits that contain the string
git log -Sdelete.me			# filter only commits adding or removing codes matching string

# checkout other commits, the files and dir will be as you commited it
git checkout b7119f2 # b7119f2 is the name of previous commit, extracted from git log output
git checkout master  # go back to working on the most recent commit
git checkout -- <file> # delete all changes done to file since last commit PERMENANTLEY






###	CHANGING HISTORY	###

# To check the difference between whats in your working dir and whats in the staging area
git diff	  # shows only changes that are still unstaged
git diff --staged # check what has been staged with git add

# create a new branch
git checkout -b <branch_name>

# see all branches created
git branch
git branch -d <branch>	# delete branch

# generate a new commit that undoes all changes introduced and apply it to new branch
# ; undoes a single commit, but after reverting this one still remains in project
# history.
git revert <commit>

# reset to older version and permenantely delete everything that had proceeded it
git reset 		# reset staging area to match the most recent commit
git reset HEAD <file>	# unstage a file (which has been git add) to same state as last commit
git reset --hard 	# as above and all uncommitted changes are gone forever
git reset --hard <commit>	# as above but one can choose a specific commit

# removing files that are untrackable in directory
git clean -n # dry run; shows what will be removed after executing
git clean -df # do it for realsies (both directories & files


	git reset --hard && git clean -df # great when everything has went wrong










###	GIT REMOTE	###

# create an alias/shorthand of a remote handle (github repo)
git remote add testing https://github.com/ultraDross/test # testing = github link
git remote rename testing development	# alter shorthand name
git remote rm development	# remove the remote

# show remote handles
git remote 	# shorthand displayed
git remote -v	# shorthand and link
git remote show testing   # see more info about testing remote
git remote add <shortname> <remote_URL> # add a new remote to local project

# download data from the remote project you dont have yet and store as a new branch,
# which could be merged.
# DIFFERENCE between clone is that fetch fetches any new work pushed since it was cloned
# or last fetched from.
git fetch testing	# branch: testing/master

# automaticatically fetch and merge a remote branch ino current branch
git pull testing

# upload your master branch to your github page
git push testing master	# only works if your work is derived from the most recent clone






###	TAGGING	###

# created a tag
git tag -a v0.1 -m "this is version 0.1"	# annotated tag contains info
git tag -a v0.1 <commit.checksum>		# annotate tag of earlier commit
git tag v0.1-light				# lightwieght tag containing less info

# list all tags, give all info of a tag
git tag -l					# list tags
gits show <tag>					# detail info about commit associated with the tag

# delete a tag
git tag -d v0.1		# delete tag v0.1

# the push command has to done seperately for the tags
git push testing --tags

# checkout a specific tagged commit on a seperate branch
git checkout -b <branch_name> <tag_name>







### BRANCHING ####
# HEAD essentially tells one which Branch they are currently working on

# reasons to use a brach: perhaps you have discovered a bug. One can start a new
# branch to fix the bug and test it. Once it is fixed, the new branch can be merged
# with the master branch

# create a new branch
git branch hotfix		# branch named hotfix produced

# move to another branch
git checkout -b master		# will only work if everything is commited in current branch

# list branches
git branch
git branch -v 		# lists more details
git branch --merged	# list branches which have been merged into the branch youre on
git branch --no-merged


# merge the tested hotfix with master branch
## if the development history has diverged from a much older point then git seems to create
## a new snapshot resulting from the three way merge; two branches snapshot and its common
## ancestor snapshot
git checkout master; git merge hotfix

# delete a branch
git branch -d hotfix

# Merge Conflicts; if you changed the same part of the same file differently in the two
# branches youre merging together, Git wont be able to merge them cleanly. It will notify
# you of the conflic and pause the process untill it is resolved (check the exact conflict
# using git status). git adds conflict resolution markers to the files so you can open
# the file and manually alter them:
#	anything marked with HEAD in the conflict file is direved from the branch which
#	you ran the merge command, everything below the ====== is from the other branch.
# you must choos one side or the other, or alter the file completely, then add/commit; thus
# resolving the merge conflict.

# graphical interface tool to fix merge conflicts
git mergetool	# needs to be confiquired or something

# suggest having master branch as 'stable' and 'developemnt' branch for new-features testing
# then when it is stable just merge with the master branch. Another branch, 'topic', is useful
# as a shrtlived branch that you use for a single feature




### REMOTE BRANCHES	###
# get remote references (branches, tags etc currently in guthub repo)
git ls-remote testing

# to sync the master (local) & origin/master (remote) branch; looks up the server 'origin'
# and fetches any data from it you don't have and updates your local database.
git fetch origin # creates a NON-editable branch, to do this; git merge origin/master

# make a local branch based off the remote branch origin/testing (called a tracking branch)
git checkout testing_remote origin/testing

# list tracking branche set up
git branch -vv

# if to set a local branch as a remote branch you just pulled down
git branch -u origin/testing

# deleting remote branches
git push origin --delete testing
