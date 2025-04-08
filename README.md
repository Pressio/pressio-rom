
# Overview

Pressio is an open-source computational framework aimed at advancing
the field of reduced-order models (ROMs) for dynamical systems in science and engineering.
We employ generic programming, and target shared
and distributed-memory applications using arbitrary data-types and diverse programming models.

We believe that the key to develop such a capability from the ground up
is to properly identify the building blocks and modularize them accordingly.
This is the approach adopted here: the library is composed of
several components that can be used independently and in a self-contained fashion,
but as a whole constitute the stack foundation of the rom component.
One of the main outcome of this is that regardless of your interest in ROMs,
you might find useful some of the components of the library.


Click below to checkout the documentation:

<a href="https://pressio.github.io/pressio-rom" target="_blank">
    <img src='./logos/logo-readme.svg' width='97%'>
</a>

## Development, versioning and backward compatibility

Until we reach a stable 1.0, please be patient and do not be surprised if the API
and functionalities somewhat rapidly change from one 0.x version to the next,
thus affecting backward compability.
Some components of pressio are more mature and stable than others,
but until we can claim the *same stability* level for all,
please keep this mind.

## Questions?
Find us on Slack: https://pressioteam.slack.com and/or
open an issue on [github](https://github.com/Pressio/pressio).


## License and Citation

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

The full license is available [here](https://pressio.github.io/various/license/).

At some point we plan to publish this, for now we have an arXiv preprint at: https://arxiv.org/abs/2003.07798.




<!-- [![Codecove](https://codecov.io/gh/Pressio/pressio/branch/master/graphs/badge.svg?precision=2)](https://codecov.io/gh/Pressio/pressio/branch/master) -->
