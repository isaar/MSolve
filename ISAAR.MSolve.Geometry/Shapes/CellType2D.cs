using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Geometry.Shapes
{
    /// <summary>
    /// Defines the shape of a cell only. Since there are no dependencies, it can also be used to map corresponding cell/element 
    /// types from one module to another.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public enum CellType2D
    {
        // 3 ---- 2							
        // |      |                            
        // |      |                            
        // 0 ---- 1
        Quad4,

        // 3 -- 6 -- 2
        // |         |
        // 7         5                            
        // |         |                            
        // 0 -- 4 -- 1                                                      
        Quad8,

        // 3 -- 6 -- 2
        // |    |    |
        // 7 -- 8 -- 5                            
        // |    |    |                            
        // 0 -- 4 -- 1                            
        Quad9,

        //    2
        //   /  \
        //  /    \                           
        // 0 ---  1                           
        Tri3,

        //     2
        //    /  \
        //   5    4
        //  /       \                            
        // 0 -- 3 -- 1                                                  
        Tri6
    }
}
