# C++ Flat Containers Library

Fast+efficient associative containers using sorted arrays, with an interface based on standard containers.

This was part of a C++ standard proposal and I recommend that you read it for more details: http://pubby.github.io/proposal.html

#### Container Adaptors
- `fc::flat_map`
- `fc::flat_multimap`
- `fc::flat_set`
- `fc::flat_multiset`

#### Class Aliases
- `fc::vector_map`
- `fc::vector_multimap`
- `fc::vector_set`
- `fc::vector_multiset`

#### New Member Functions
- `has`
   - `map.has(key)` returns a pointer to key's mapped value if it exists, otherwise returns null.
-  `insert` (delayed sort)
    - `map.insert(first, last, fc::delay_sort)` performs insertion with a delayed sort optimization.
- Constructors (delayed sort overload)
    - Constructs flat container using a delayed sort optimization.
- Constructors (container construct overload)
    - Constructs flat container by forwarding arguments to the underlying container member.

#### What's an adaptor?

Container adaptors allow you to use any vector-like sequence container for the underlying storage.
`std::vector` is the natural choice, but classes like `boost::small_vector` and `boost::static_vector` are useful too. Because `std::vector` is the most common in usage, aliases are provided for convenience.

For basic use, just use the `std::vector` aliases.

#### The public Members: `container` and `underlying`

The public member `container` allows access to the underlying storage container. Note that it is the user's responsibility to keep the container in a sorted, unique state.

Likewise, the public member of iterators: `underlying`, allows access to the underlying storage container's iterators.

*Example: Delayed sort optimization using `.container`*

    std::flat_multiset<int> set;
    while(cpu_temperature() < 100)
        set.container.push_back(cpu_temperature());
    std::sort(set.begin(), set.end());

#### Const Iteration by Default

For safety reasons, flat container iterators are const by default. To bypass this safety and get non-const iterators, one can either iterate `.container` or take the `.underlying` of the iterator.

*Example: Modify values in a way that preserves sortedness*

    for(auto& v : set.container)
        v *= 2;

*Example: Same thing but with `.underlying`*

    for(auto it = v.begin(); it != v.end(); ++it)
        (*v.underlying) *= 2;

#### Helper Types

The directory `include_extra` contains convenience typedefs for use with Boost.Container.
