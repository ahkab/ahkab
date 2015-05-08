Jekyll Clean
============

* Get it from [github](https://github.com/scotte/jekyll-clean).
* See the [live demo](https://scotte.github.io/jekyll-clean).
* See it [in action on my own blog](https://scotte.github.io).

A simple and clean Jekyll theme using [bootstrap](http://getbootstrap.com)
(not to be confused with jekyll-bootstrap) that's easy to modify and very
modular in component and element reuse.

It uses Disqus for comments and includes Google Analytics support. Both of
these features are disabled by default and can be enabled via \_config.yml. You
can also rip this code out of the templates if you like (footer.html and post.html).
The beauty of Jekyll - keep things clean... Jekyll Clean!

The theme works well on mobile phones, using a collapsable nav bar and hiding the
sidebar. The links pane in the sidebar is available on mobile through the nav menu,
and you can do the same thing for any other sections added to the sidebar.

Don't forget to occassionally merge against my upstream repository so you can get
the latest changes. Pull requests are encouraged and accepted!

Installation
============

If you don't have a blog already on github, start by cloning this repository.
Best to do that directly on github and then clone that down to your computer.

If you already do have a blog, You can certainly apply this theme to your existing
blog in place, but then you won't be able to merge as the theme changes. If you
re-apply your blog history on top of this theme's **gh-pages** branch, it's then
easy to update to the latest version of the theme. You also don't want to have to
deal with resolving old conflicts from your existing history, so you may wish to to
push your existing master off to a new branch so you have the old history and start
a new branch with this as the start, merging in your \_posts and other assets (after
git rm'ing the current \_posts.

Not ideal, but you have to make a choice - either apply it manually or base your
blog off this theme's branch. Either way it will work, and both have their own
pros and cons.

You can setup an upstream tracking repository like so:

```
$ git remote add upstream git@github.com:scotte/jekyll-clean.git
```

And now when you wish to merge your own branch onto the latest version of the
theme, simply do:

```
$ git fetch upstream
$ git merge upstream/gh-pages
```

Of course you will have to resolve conflicts for \_config.yml, \_includes/links-list.html,
and \_posts, and so on, but in practice this is pretty simple.

This is how I maintain my own blog which is based on this theme. The old history is
sitting in an **old-master** branch that I can refer to when I need to.

License
=======

The content of this theme is distributed and licensed under a
![License Badge](/images/cc_by_88x31.png)
[Creative Commons Attribution 4.0 License](https://creativecommons.org/licenses/by/4.0/legalcode)

    This license lets others distribute, remix, tweak, and build upon your work,
    even commercially, as long as they credit you for the original creation. This
    is the most accommodating of licenses offered. Recommended for maximum
    dissemination and use of licensed materials.

In other words: you can do anything you want with this theme on any site, just please
provide a link to [the original theme on github](https://github.com/scotte/jekyll-clean)
so I get credit for the original design. Beyond that, have at it!

This theme includes the following files which are the properties of their
respective owners:

* js/bootstrap.min.js - [bootstrap](http://getbootstrap.com)
* css/bootstrap.min.css - [bootstrap](http://getbootstrap.com)
* js/jquery.min.js - [jquery](https://jquery.com)
* images/cc_by_88x31.png - [creative commons](https://creativecommons.org)
