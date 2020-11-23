## Setting up git lfs

[git lfs](https://git-lfs.github.com/) is used to support large files - usually data files which are too large to track . Here we recommend that you let git lfs manage all files inside received/. To do this, you only need run a couple of git lfs commands, as shown in [here](https://git-lfs.github.com/).

```
> git lfs install
> git lfs track 'received/**'
```

To see the file which are tracked by git lfs, use

```
> git lfs ls-files
```