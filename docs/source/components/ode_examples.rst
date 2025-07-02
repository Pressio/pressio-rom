example 1
=========

.. literalinclude:: ../../../tests/functional_small/ode_steppers/doc_main1.cc
   :language: cpp
   :lines: 48-


example 2
=========

Integrate the `original Lorenz system <https://en.wikipedia.org/wiki/Lorenz_system>`__ via explicit time integration.

.. literalinclude:: ../../../tests/functional_small/ode_steppers/doc_main2.cc
   :language: cpp
   :lines: 48-


Visualize
---------

You can use the data collected to plot the results via the following script: 


.. code-block:: py

   import numpy as np
   import matplotlib.pyplot as plt
   from mpl_toolkits.mplot3d import Axes3D

   if __name__== "__main__":
       numVars = 3
       D = np.fromfile("./state_snapshots.bin")
       numTimeSteps = int(np.size(D)/numVars)
       print(numTimeSteps)
       #D  = np.reshape(data, (numTimeSteps, numVars))

       fig = plt.figure()
       ax = fig.add_subplot(1, 1, 1, projection='3d')
       ax.plot(D[0::3], D[1::3], D[2::3])
       fig.savefig("plot.png", format="png", bbox_inches='tight', dpi=300)
       plt.show()


.. image:: ../_static/ode_ex1.png
  :width: 50 %
  :align: center
  :alt: Lorenz


example 3
=========

Integrate the `original Lorenz system <https://en.wikipedia.org/wiki/Lorenz_system>`__ using implicit time integration.

.. literalinclude:: ../../../tests/functional_small/ode_steppers/doc_main3.cc
   :language: cpp
   :lines: 48-

Visualize
---------

You can visualize the results using the same Python script above.
