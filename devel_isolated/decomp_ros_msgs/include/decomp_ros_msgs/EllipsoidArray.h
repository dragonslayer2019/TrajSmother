// Generated by gencpp from file decomp_ros_msgs/EllipsoidArray.msg
// DO NOT EDIT!


#ifndef DECOMP_ROS_MSGS_MESSAGE_ELLIPSOIDARRAY_H
#define DECOMP_ROS_MSGS_MESSAGE_ELLIPSOIDARRAY_H


#include <string>
#include <vector>
#include <memory>

#include <ros/types.h>
#include <ros/serialization.h>
#include <ros/builtin_message_traits.h>
#include <ros/message_operations.h>

#include <std_msgs/Header.h>
#include <decomp_ros_msgs/Ellipsoid.h>

namespace decomp_ros_msgs
{
template <class ContainerAllocator>
struct EllipsoidArray_
{
  typedef EllipsoidArray_<ContainerAllocator> Type;

  EllipsoidArray_()
    : header()
    , ellipsoids()  {
    }
  EllipsoidArray_(const ContainerAllocator& _alloc)
    : header(_alloc)
    , ellipsoids(_alloc)  {
  (void)_alloc;
    }



   typedef  ::std_msgs::Header_<ContainerAllocator>  _header_type;
  _header_type header;

   typedef std::vector< ::decomp_ros_msgs::Ellipsoid_<ContainerAllocator> , typename std::allocator_traits<ContainerAllocator>::template rebind_alloc< ::decomp_ros_msgs::Ellipsoid_<ContainerAllocator> >> _ellipsoids_type;
  _ellipsoids_type ellipsoids;





  typedef boost::shared_ptr< ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator> > Ptr;
  typedef boost::shared_ptr< ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator> const> ConstPtr;

}; // struct EllipsoidArray_

typedef ::decomp_ros_msgs::EllipsoidArray_<std::allocator<void> > EllipsoidArray;

typedef boost::shared_ptr< ::decomp_ros_msgs::EllipsoidArray > EllipsoidArrayPtr;
typedef boost::shared_ptr< ::decomp_ros_msgs::EllipsoidArray const> EllipsoidArrayConstPtr;

// constants requiring out of line definition



template<typename ContainerAllocator>
std::ostream& operator<<(std::ostream& s, const ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator> & v)
{
ros::message_operations::Printer< ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator> >::stream(s, "", v);
return s;
}


template<typename ContainerAllocator1, typename ContainerAllocator2>
bool operator==(const ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator1> & lhs, const ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator2> & rhs)
{
  return lhs.header == rhs.header &&
    lhs.ellipsoids == rhs.ellipsoids;
}

template<typename ContainerAllocator1, typename ContainerAllocator2>
bool operator!=(const ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator1> & lhs, const ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator2> & rhs)
{
  return !(lhs == rhs);
}


} // namespace decomp_ros_msgs

namespace ros
{
namespace message_traits
{





template <class ContainerAllocator>
struct IsFixedSize< ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator> >
  : FalseType
  { };

template <class ContainerAllocator>
struct IsFixedSize< ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator> const>
  : FalseType
  { };

template <class ContainerAllocator>
struct IsMessage< ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator> >
  : TrueType
  { };

template <class ContainerAllocator>
struct IsMessage< ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator> const>
  : TrueType
  { };

template <class ContainerAllocator>
struct HasHeader< ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator> >
  : TrueType
  { };

template <class ContainerAllocator>
struct HasHeader< ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator> const>
  : TrueType
  { };


template<class ContainerAllocator>
struct MD5Sum< ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator> >
{
  static const char* value()
  {
    return "e2c31e58d2b4b09679be4a3c12fffb19";
  }

  static const char* value(const ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator>&) { return value(); }
  static const uint64_t static_value1 = 0xe2c31e58d2b4b096ULL;
  static const uint64_t static_value2 = 0x79be4a3c12fffb19ULL;
};

template<class ContainerAllocator>
struct DataType< ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator> >
{
  static const char* value()
  {
    return "decomp_ros_msgs/EllipsoidArray";
  }

  static const char* value(const ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator>&) { return value(); }
};

template<class ContainerAllocator>
struct Definition< ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator> >
{
  static const char* value()
  {
    return "Header header\n"
"Ellipsoid[] ellipsoids\n"
"\n"
"================================================================================\n"
"MSG: std_msgs/Header\n"
"# Standard metadata for higher-level stamped data types.\n"
"# This is generally used to communicate timestamped data \n"
"# in a particular coordinate frame.\n"
"# \n"
"# sequence ID: consecutively increasing ID \n"
"uint32 seq\n"
"#Two-integer timestamp that is expressed as:\n"
"# * stamp.sec: seconds (stamp_secs) since epoch (in Python the variable is called 'secs')\n"
"# * stamp.nsec: nanoseconds since stamp_secs (in Python the variable is called 'nsecs')\n"
"# time-handling sugar is provided by the client library\n"
"time stamp\n"
"#Frame this data is associated with\n"
"string frame_id\n"
"\n"
"================================================================================\n"
"MSG: decomp_ros_msgs/Ellipsoid\n"
"float64[3] d\n"
"float64[9] E\n"
;
  }

  static const char* value(const ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator>&) { return value(); }
};

} // namespace message_traits
} // namespace ros

namespace ros
{
namespace serialization
{

  template<class ContainerAllocator> struct Serializer< ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator> >
  {
    template<typename Stream, typename T> inline static void allInOne(Stream& stream, T m)
    {
      stream.next(m.header);
      stream.next(m.ellipsoids);
    }

    ROS_DECLARE_ALLINONE_SERIALIZER
  }; // struct EllipsoidArray_

} // namespace serialization
} // namespace ros

namespace ros
{
namespace message_operations
{

template<class ContainerAllocator>
struct Printer< ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator> >
{
  template<typename Stream> static void stream(Stream& s, const std::string& indent, const ::decomp_ros_msgs::EllipsoidArray_<ContainerAllocator>& v)
  {
    s << indent << "header: ";
    s << std::endl;
    Printer< ::std_msgs::Header_<ContainerAllocator> >::stream(s, indent + "  ", v.header);
    s << indent << "ellipsoids[]" << std::endl;
    for (size_t i = 0; i < v.ellipsoids.size(); ++i)
    {
      s << indent << "  ellipsoids[" << i << "]: ";
      s << std::endl;
      s << indent;
      Printer< ::decomp_ros_msgs::Ellipsoid_<ContainerAllocator> >::stream(s, indent + "    ", v.ellipsoids[i]);
    }
  }
};

} // namespace message_operations
} // namespace ros

#endif // DECOMP_ROS_MSGS_MESSAGE_ELLIPSOIDARRAY_H