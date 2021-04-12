# Contributing

:construction: Work in Progress :construction:
## Setting up your fork

First fork the repo on github.  When you have created your fork, clone it so that you have a local copy:

```sh
git clone https://github.com/your-github-username/pyfolding.git
cd pyfolding
```

Next, you need to link it to the main repo:
```sh
git remote add upstream https://github.com/quantumjot/pyfolding.git
```
**This is important**, because it enables you to keep your fork up-to-date with everyone's changes.

## Making a new branch

If your want to contribute your new feature, make sure to make a new branch:
```sh
git checkout master -b your-branch-name
```

Use the following to determine which branch you're on:
```sh
git status
```

## Keeping your fork up to date
The following commands will update the main branch of your fork:
```sh
git checkout master
git pull upstream master
```

You can then update any branches by:
```sh
git checkout your-branch-name
git merge master
```
