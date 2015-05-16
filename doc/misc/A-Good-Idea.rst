Why `ahkab`
===========

A rather long winded explanation as to why this project exists today, you may
skip directly to :ref:`install-ahkab`.

Do circuit simulators dim your wit?
-----------------------------------

A young engineer begins to design a well-known, basic circuit block. He
knows his stuff and he's prepared: he did his homework on pen and paper
first and now he wants to take care of schematic entry, checking that the
circuit indeed works as expected, in case introducing slight adjustments as
needed.

.. image:: ../images/mixer-sb.png
   :alt: A basic circuit block
   :align: center

He draws the circuit and, without even thinking about it, he clicks *Simulate*.

The results are surprising to him, something he has not taken into account has
influenced his simulation and the numbers are slightly off. Or maybe the
schematic was drawn in a rush?

He starts fiddling with the design parameters. Make some transistors
bigger. Now, make a few others smaller. *Shouldn't this fix it?* - he
wonders. *Maybe*, but the aspect ratios that were just changed almost at
random were initially selected to reach multiple results at the same
time.

Now, the design does not meet, not one, but multiple specs. He removes parts
that are secondary for the current simulation. Then proceeds to change the
simulation itself...

What's left of a good idea
''''''''''''''''''''''''''

If a fellow design engineer were to look at his screen now, he'd see a
horrifying circuit: aspect ratios all over the place, cheap hacks to
make up for the "secondary" parts that were removed in order to simplify
the circuit, quantities that make no sense in the physical world.

After going down the road of compulsive circuit simulation, little is
left of the initial, promising design and our young engineer feels lost,
frustrated and he's metaphorically about to hitting his head against the
keyboard of his workstation.

*Where's the fun of doing circuit design this way?*

In the end
''''''''''

We have seen the above before... and that's not how our beloved
microelectronic work should be.

Circuit simulators are just a tool. An insidious one at it, as we may naively
use them instead of our gray matter, blindly trusting our models or giving in to
the temptation of manual, broken "optimization" fiddling, rather than as a
verification tool, area where they truly excel.

It doesn't help much that often neither it is straightforward to debug a
simulation, nor it is clear what exactly the simulator is doing.

Which brings us to:

What we try to do here with ``ahkab``
-------------------------------------

The ``ahkab`` circuit simulator is an experiment.

We have no expectation that our proof-of-concept, sometimes buggy, small
circuit simulation tool will be replacing the mainstream circuit simulators:
they are mainstream for good reasons and they very much deserve the praise and
money we pay. It would be foolish to think otherwise.

1. Peek under the hood
''''''''''''''''''''''

But we do still think we have our own place: what we wish to do with
``ahkab`` is to allow the user, the designer, to peek behind the veil
and see more clearly what goes on with his simulations.

For this reason, ``ahkab`` supports operations such as printing out all
equation matrices.

It is also written in a scripted, interpreted language (Python) that,
while requiring us to sacrifice raw speed, should make it relatively
easy to see what's going on behind the hood.

And the algorithms are there for you to see, inspect and, if need be,
correct: too often scientific papers about software come with no
available implementation or the source is not distributed: with
``ahkab`` all code is available under a copy-left license, allowing you
to benefit of the code, modify it and giving others the same freedom.

2. Experiment
'''''''''''''

With ``ahkab``, you are welcome to implement an algorithm you have read
about in a paper, if you are so inclined: we believe no lecture, no
matter how in-depth, will provide an insight in a circuit simulation
algorithm as deep as rolling up your sleeves and trying to code it up
yourself. (*please remember talk to us well in advance if you expect us
to include your work!*)

3. Have fun doing electronics!
''''''''''''''''''''''''''''''

All in all, we hope this little project helps you understand better what
will goes on in your circuit when you implement it, when you simulate it
and especially *we wish you have fun while doing so!*

