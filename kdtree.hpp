#pragma once

#include <type_traits>
#include <list>
#include <vector>



struct Splitter /*concept*/
{
    // split the given bounding box
    // given bounding box will then be the lower part, and the returned bounding box will be the upper part
    template<class real>
    BoundingBox split(BoundingBox<real>& bb, Points& points /* maybe it would also be usefull to give the last dimension that has been split*/);

};

/*!
 * \bried Define the plane that splits space.
 *        It holds only two variables: the dimension n the plane is perpendicular to
 *        and the anchor point of the plane, i.e. the nth coordinate of a point on the plane.
 */
template<class real>
struct HyperPlane
{
    std::size_t dimension; //<! which dimension is the plane perpendicular to
    real        anchor;    //<! nth component of coordinate of a point on the plane
};

/*!
 * \brief Implementation of Splitter concept
 * Subdivide each cell into two equally sized cells, while rotating the dimension of the splitting hyperplance
 */
struct RotatingSubdivision
{
    // we can access the points contained in node via the nodes iterator
    // it is assumed, that the node is not yet split
    // we need to get information about the splitting of the nodes parent
    template<class BoundingBox, class Node>
    void split(BoundingBox& bb, const Node& node)
    {
        // the dimension that is to be split
        std::size_t d = (node.parent().spitting_dimension() + 1) % bb.dimension();
        // the location of the hyperplane TODO use correct coordinate type
        double c = bb.lower[d] + (bb.upper[d]-bb.lower[d]) / 2.;
        // or return Seperator(d,c);
        return std::make_tuple(d,c);
    }
};

/*!
 * \brief Provide default interface for value_type to coordinate transformation
 */
struct DefaultA
{
    // deduce the real type of the coordinates, strip all possible references and consts
//     typedef typename std::remove_const<typename std::remove_reference<decltype(std::declval<const V&>().operator[](0))>::type>::type real;

    // access the i.th dimension of the coordinate of v
    template<class V>
    auto operator()(const V& v, std::size_t t) const -> decltype(v[i]) { return v[i];}

    // get the dimension of space v is located within
    // TODO: because we assert that all points live in the same space, this should always evaluate to the same value
    // and is hence somewhat redundant, maybe it would be better to give the dimension to the tree uppon construction?!?
    template<class V>
    std::size_t size(const V& v)               const { return v.size();}
};


/*!
 * \tparam V the value type that is to be sorted into the tree
 * \tparam NV data that is associated with nodes not leafes
 * \tparam S the splitting rule for subdivision of nodes
 * \tparam A an object that is capable to obtain coordinates from V and the dimension of V
 * #\tparam VC the container used to store value types (defaults to std::list)
 * #\tparam NC the container used to store the nodes (defaults to std:list)
 *
 */
template<class V, class NV, class S = RotatingSubdivision, class A = DefaultA/*, class VC, class NC*/>
class kdtree
{
public:
    typedef typename Traits::real_type real_type;

    typedef typename Traits::Splitter Splitter;

    typedef typename Traits::Coordinate Coordinate;

    typedef typename Traits::value_type value_type;

    typedef std::list<Node> NodeContainer;

    typedef std::list<value_type> DataContainer;

    typedef typename NodeContainer::iterator node_iterator;

    typedef typename Container<value_type>::iterator data_iterator;

    struct BoundingBox
    {
        real_type* _lower; // array holding the coordinates of the lower corner
        real_type* _upper; // array holding the coordinates of the lower corner
    };

    // a node is considered a leightweight, iterator like thing
    // it does not store user data, and gets constructed by the tree itself.
    // it provides functionality to iterate over all contained data and all contained subnodes
    // it only holds the splitting dimension and position, no real boundingbox
    // if the user requires more functionality (e.g. attached data to non leaf nodes),
    // he can provide a NodeExtension class with the traits.
    // the NodeBase class constructor will get the boundingbox of the node and the iterator range
    // for the data that will be contained within this node
    // a possibility would also be to store a user defined node specific data type, but with inheritance,
    // we allow the compiler to use empty base class optimization, resulting in zero overhead if the Node specific
    // data type is independent.
    // maybe we can clarify this by using private inheritance and a node data getter that returns a cast of this to
    // user data type (see specific commented lines)
    // struct Node : private Nodebase
    
    // provide access user node data
    struct UserNodeBase {
        // access data
        const NV& data() const { }
    };
    
    // do NOT provide access to user node data
    struct NoUserNodeBase{};
    
    // choose base class type according to whether there is a user Node data or not
    using NodeBase = std::conditional<std::is_same<NV,EmptyNode>::value>,NoUserNodeBase,UserNodeBase::type;
    
    // inherit from the correct Node Base, i.e. add a data member if NV is not the default value,
    // otherwise inherit from NoUserNodeBase and use empty base class optimization
    struct Node : public NodeBase, public HyperPlane<real_type>
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
        // TODO decide: access splitter information by base type? or provide direct access?
        // std::size_t splitting_dimension() const;
        // real_type   splitting_position() const;
        const HyperPlane<real_type>& splitter() const
        {
            return static_cast<const HyperPlane<real_type>&>(*this);
        }

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

    DataContainer<value_type> _data;  // store all data in a linear container. likely a std::vector is usefull, because we need random access and swaping items
    NodeContainer<Node>       _nodes; // store all nodes in a linear container. likely a std::list is usefull, because building the tree needs lots of insertions at arbitrary position
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
    
    // moving a point will probably require a complete rebuild of the tree because
    // we do not now if the old splitting is still consistent with the moved point(s)
    // maybe we should not allow that at all...
    void move_data(...);

    // access the tree for searching and other stuff
    Node& root();
    const Node& root() const;

    // iterate overall nodes, the order is given by the internal order of the node container,
    // i.e. first node will be the root, then the left subtree will be iterated (recursively)
    // last the right subtree will be iterated
    // TODO: it would be nice to have a range returned, but that would add dependency on boost :-/
    //       maybe one can enable this feature with #ifdef if boost defines are available?!?
    // boost::iterator_range<> nodes() const;
    node_iterator node_begin() const;
    node_iterator node_end() const;

    // iterate over data the order is again defined by the tree structure,
    // the left most data elements are first, then the leafes are iterated towards the rightmost leaf of the tree
    // TODO: it would be nice to have a range returned, but that would add dependency on boost :-/
    //       maybe one can enable this feature with #ifdef if boost defines are available?!?
    // boost::iterator_range<> data() const;
    // TODO think about constness, obviously sometimes one will want to modify some of the contained data.
    //      WARNING: modifying contained data, i.e. the datas location in space will invalidate the tree!
    data_iterator data_begin() const;
    data_iterator data_end() const;
};

