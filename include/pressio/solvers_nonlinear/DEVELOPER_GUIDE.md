The design of nonlinear solvers is organized into three main components:
(1) a public API, (2) an internal representation of the coordination "loop" to handle the steps executing the solve, and (3) a collection of internal free functions.
This layered structure aims to separate concerns, improve extensibility, and make the solver components modular and reusable.

The public API for nonlinear solvers in `pressio-rom` is defined across several headers in the `include/pressio/solvers_nonlinear` directory.
It provides a set of free functions that users call to create solver objects. These functions are the only interface that users are expected to interact with.
Each function corresponds to a specific solver type, such as Newton-Raphson or Gauss-Newton. The user calls one of these functions with the required arguments. The function returns a solver object whose type is intentionally hidden from the user. The user does not need to know the internal type of the solver object; they only need to know that the object exposes specific methods, such as `solve(system, state)`.
Internally, when a create function is called, the code first constructs the appropriate registry for the selected solver type. The registry is configured to store the temporary data needed by the solver, such as residual vectors and Jacobians. Next, the code composes this registry with the solver’s iteration loop implementation. Finally, the full solver object is instantiated and returned to the user.

The registry is an internal component used by nonlinear solvers in `pressio-rom` to store and manage data that is shared across different parts of the solver. It is implemented in the file `impl/registries.hpp` and plays a central role in decoupling the components of the solver loop.
The main purpose of the registry is to hold data that is needed during the solve process. Examples of such data include the residual vector, Jacobian matrix, update direction, and any other temporary quantities used during the nonlinear iterations. By storing these in a central location, the solver avoids having to pass them explicitly between all the internal functions that participate in the solver loop.
Each solver type has its own registry. The registry internally stores data using member variables, typically with types deduced from the system and state. The exact contents of the registry depend on the solver type. For example, a Newton-Raphson solver registry might store a residual vector, a Jacobian matrix, and an update vector, while a Levenberg-Marquardt solver registry might also include a damping factor or additional vectors.
The registry is created internally when the user constructs a solver through the public API.
Once created, the registry is passed into the solver loop and used throughout the nonlinear iterations.
The interal "solve" loop/coordinator takes care of calling internal free functions accessing and modifying data in the registry as needed.
This design offers several advantages. It reduces the number of arguments passed around in internal functions, which improves readability and flexibility.
It also allows different components to operate independently, as long as they agree on what is stored in the registry and under what name or interface.
Importantly, the registry is not exposed to users. It is created, owned, and used entirely within the solver infrastructure.
Users constructing solvers and calling `solve` are not aware of the registry’s existence, and they are not expected to interact with it.
This registry approach is just an implementation detail.

How does this registry thing work?

To keep things organized and type-safe, the registry uses a system of tags.
A tag is just a small type that acts like a label for a specific kind of data.
Instead of using strings or variable names, each kind of data is associated with a unique tag. For example, there might be a `ResidualTag` for the residual vector and a `JacobianTag` for the Jacobian matrix. These are small struct types defined specifically for labeling purposes.
When the registry is created, it constructs all the required data objects and assigns each of them to its corresponding tag. This means that at the time the registry object is initialized, all necessary internal data (like residual, Jacobian, etc.) is also allocated and stored, specifically suited for the problem being solved. This construction process is specific to each solver type and ensures that the registry has everything it needs before the solver loop begins.
The data is stored as member variables inside the registry, and each tag maps to a specific one of these variables. The type of each value is fixed and known based on the tag.
One additional benefit of using tags is that they make it easy to understand what data each solver uses.
Each solver has its own registry definition, and the tags used in the registry are clearly listed in `registries.hpp`. Because these tags have expressive names—like `ResidualTag`, `JacobianTag`, or `UpdateDirectionTag`—a developer can simply look at the list of tags in the registry to get a quick and complete view of all the internal data that the solver relies on. This is helpful for both understanding and debugging the solver’s behavior, and also when extending the implementation or adding a new solver type.
To access data from the registry, functions use the tag. For example, if a function needs the residual vector, it can ask the registry to return the object associated with the `ResidualTag`. The registry looks up the correct member and returns a reference to it.
This system allows functions to access shared data without needing to know how the registry is structured. Each function just needs to know which tag to use.
The advantage of this approach is that it avoids using strings or dynamic lookups. It is checked at compile time, so mistakes like asking for the wrong type of object are caught early. It also helps keep the code modular, because different parts of the solver can work independently as long as they agree on which tags to use.

The design of using a registry, an implementation loop, and free functions acting on tagged data is particularly well suited for nonlinear solvers because it promotes a clean separation of concerns, flexibility, and performance-aware design. Nonlinear solvers are often composed of small, distinct operations—such as computing a residual, forming a Jacobian, solving a linear system, updating the state, and checking convergence. These operations have clear boundaries and are typically stateless. Modeling them as free functions rather than member functions of large classes allows them to be implemented, tested, reused, and composed more easily.
The registry approach helps manage shared data between these operations without tightly coupling them.
Each component only needs to know what data it uses (via tags), not how or where it is stored.
This makes the code modular and reduces the need to pass large numbers of arguments around.
It also allows the solver loop to be written generically, without knowing in advance what exact data is involved.
The solver loop acts as the coordinator. It knows the sequence of operations but delegates the actual computations to free functions and uses the registry to access the shared data.
This design is especially helpful when implementing multiple solver types (e.g., Newton-Raphson, Gauss-Newton, Levenberg-Marquardt), because the loop can be customized or reused with minimal duplication.


The file `functions.hpp` defines a set of internal helper functions that carry out the basic operations needed during a nonlinear solve. These functions are not solver-specific—they are used by different solver types and are called from within the solver loop. Each function is responsible for a focused task, such as computing a residual, forming a Jacobian, updating the state.
These functions are designed to operate on data stored in the registry, using tag-based access, and they rely on the system and state provided during the solve.
By organizing the logic into separate free functions, the code stays modular, easier to test, and easier to extend.
For example, here is a high-level description of the roles of some key functions defined in this file:
* **compute\_residual**: Computes the nonlinear residual vector by evaluating the user-provided system at the current state and storing the result in the registry.
* **compute\_jacobian**: Computes the Jacobian matrix of the system at the current state and stores it in the registry.
* **compute\_residual\_and\_jacobian**: A combined version that evaluates both the residual and the Jacobian in a single call, if the system supports doing so efficiently.
These functions form the building blocks of the solver loop. Each iteration of the solver calls some or all of them, depending on the solver strategy.
Because each function performs a single task, the design makes it easier to adapt the loop for different solvers or plug in custom logic if needed.
They work together with the registry and the solver loop to drive the nonlinear solve.
