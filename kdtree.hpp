#pragma once

#include <list>
#include <vector>


struct Splitter /*concept*/
{
    // split the given bounding box
    // given bounding box will then be the lower part, and the returned bounding box will be the upper part
    BoundingBox split(BoundingBox& bb, Points& points /* maybe it would also be usefull to give the last dimension that has been split*/);

};

struct kdtreetraits /*concept*/
{

    // the type of data that is stored within the tree
    typedef void value_type;

    // the (scalar) numeric type used for coordinates that define the location of data in space
    typedef void real_type;

    // the type of the splitter to be used for splitting cells
    typedef void Splitter;

    // a functor that is capable to receive the Nth coordinate associated with a value_type
    // i.e. its singature is supposed to be one of the following:
    //      const real_type& operator()(const value_type& item, std::size_t dim) const
    // or.
    //      real_type operator()(const value_type& item, std::size_t dim) const
    // alternatively, because the dimensionality of the space also needs to be defined
    //      iterator_range operator()(const value_type& item) const
    // then, one can use
    //      std::distance(range.begin(), range.end())
    // to get dimensionality of the points and
    //      *(range + n)
    // to access nth coordinate
    typedef void Coordinate;

    // the linear storage type to be used for the data within the tree
    // WHY do we need Container as template parameter to the algorithm?
    // Answer, some value_types (e.g. Eigen's static size matrices) need aligned allocation.
    // Hence special care for storing those items is required.
    typedef void Container; //TODO: should this be a template? because in the end we probably use Container<value_type> or container<Node<value_type>>

    // generate a splitter object
    Splitter splitter_object();

    // generate coordinate getter
    Coordinate coordinate_object();
};

// a default model that implements the traits concept
class DefaultkdtreeTraits
{
};


// TODO: do we want to combine all template parameters within the traits class, or provide a long list of seperate template parameters
template<class Traits>
class kdtree
{
public:
    typedef typename Traits::real_type real_type;

    typedef typename Traits::Container Container;

    typedef typename Traits::Splitter Splitter;

    typedef typename Traits::Coordinate Coordinate;

    typedef typename Traits::value_type value_type;

    typedef std::list<Node> NodeContainer;

    typedef typename NodeContainer::iterator node_iterator;

    typedef typename Container<value_type>::iterator data_iterator;

    typedef typename Traits::NodeBase NodeBase; // TODO: maybe we can user SFINAE to fallback to the default implementation if the user does not define a NodeBase within his traits class?

    struct BoundingBox {};

    // a node is considered a leightweight, iterator like thing
    // it does not store user data, and gets constructed by the tree itself.
    // it provides functionality to iterate over all contained data and all contained subnodes
    // it only holds the splitting dimension and position, no real boundingbox
    // if the user requires more functionality (e.g. attached data to non leaf nodes),
    // he can provide a NodeExtension class with the traits.
    // the NodeBase class constructor will get the boundingbox of the node and the iterator range
    // for the data that will be contained within this node
    struct Node : public NodeBase
    {
        // defines the iterator range of data that is contained within this node an all of its children
        data_iterator _data_begin;
        data_iterator _data_end;

        // defines the iterator range of for nodes that is contained within this node an all of its children
        node_iterator _node_begin;
        node_iterator _node_end;

        // track children, if this is a leaf node, dereferencation is invalid
        node_iterator _lower;
        node_iterator _upper;

        Node(const BoundingBox& bb, data_iterator db, data_iterator de) :
            NodeBase(bb,db,de),
            _data_begin(db),
            _data_end(de)
            {}

        // each node splits space with a hyperplance that is perpendicular to its splitting dimension
        // hence, the hyperplane is completly defined by the index of the dimension that is split
        // and the coordinate of the plance wrt to this index.
        std::size_t splitting_dimension() const;
        real_type   splitting_position() const;

         // iterate overall children
        boost::iterator_range<> nodes() const;

        // iterate over data contained in children
        boost::iterator_range<> data() const;

        bool is_leaf();
    };

protected:

    // all user data and all nodes get stored in the tree itself. this allows putting it into a "dense package" within the memory
    // building a tree requires the data in the containers to be ordered in a specific state:
    // the root is supposed to be the first entry of the node container. children can then
    //   - be sorted depth first, the node will be succedded with the nodes of the lower subtree, then upper subtree (recursively)
    //       + easy iteration over subtrees, because they are consecutive in memory
    //       - building the tree requires shifting the subtrees in the container, but we can use e.g. std::list
    //   - be somehow different
    //       - i can not see any advantage of any other order
    //
    // the data should then be ordered in the same fashion, then each node can keep iterator ranges to data, that is associated to this node or all its children
    // one could also move this functionality into a traits class (e.g. treebuilder or similar) to allow kdtrees that user a completely different architecture
    // (e.g. for growing trees that do not require rebuild)

    Container<value_type> _data;  // store all data in a linear container. likely a std::vector is usefull, because we need random access and swaping items
    Container<Node>       _nodes; // store all nodes in a linear container. likely a std::list is usefull, because building the tree needs lots of insertions at arbitrary position
    Splitter  _splitter;
    Coordinate _coordinate;   // instance of coordinate getter
    std::size_t _bucket_size; // how many data elements are allowed per leaf, TODO: put it into the traits class? then it can also be compile time constant

public:

    kdtree(const Traits& traits = Traits()) :
        _data(),
        _nodes(),
        _splitter(traits.splitter_object()),
        _coordinate(traits.coordinate_object())
    {
    }

    // insert data, will invalidate the tree -> rebuild required
    void insert(value_type v);

    template<class input_iterator>
    void insert(input_iterator begin, input_iterator end);

    // generate internal tree structure, i.e. sort the containers and initialize nodes
    void build();

    // get the box that contains all data
    BoundingBox bounding_box() const;

    // access the tree for searching and other stuff
    Node& root();

    // iterate overall nodes, the order is given by the internal order of the node container,
    // i.e. first node will be the root, then the left subtree will be iterated (recursively)
    // last the right subtree will be iterated
    boost::iterator_range<> nodes() const;

    // iterate over data the order is again defined by the tree structure,
    // the left most data elements are first, then the leafes are iterated towards the rightmost leaf of the tree
    boost::iterator_range<> data() const;


};

