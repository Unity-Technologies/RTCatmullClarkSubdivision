using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

public static class CatmullClarkTessellation
{
	// face attributes
	public struct cct_Face
	{
		public int3 halfEdgeIDs;
	};

	// diamond parent
	public struct cct_DiamondParent
	{
		public ConcurrentBinaryTree.Node baseNode, topNode;
	};



	/*******************************************************************************
	 * FindMSB -- Returns the position of the most significant bit
	 *
	 */
	public static int cct__FindMSB(uint x)
	{
		int msb = 0;

		while (x > 1u) {
			++msb;
			x>>= 1;
		}

		return msb;
	}


	/*******************************************************************************
	 * NextPowerOfTwo -- Returns the closest power of two
	 *
	 */
	public static int cct__NextPowerOfTwo(int x)
	{
		x--;
		x|= x >> 1;
		x|= x >> 2;
		x|= x >> 4;
		x|= x >> 8;
		x|= x >> 16;
		x++;

		return x;
	}


	/*******************************************************************************
	 * BaseFaceCount -- returns the number of faces at root
	 *
	 */
	public static int cct__BaseFaceCount(CatmullClark.cc_Subd subd)
	{
		return CatmullClark.ccm_HalfedgeCount(subd.cage);
	}

	public static int cct__BaseFaceCount(CatmullClarkGPU.cc_SubdGPU subd)
	{
		return CatmullClarkGPU.ccm_HalfedgeCount(subd.cage);
	}


	/*******************************************************************************
	 * MinCbtDepth -- returns the minimum depth required for the CBT to work
	 *
	 */
	public static int cct__MinCbtDepth(CatmullClark.cc_Subd subd)
	{
		int baseFaceCount = cct__BaseFaceCount(subd);

		return cct__FindMSB((uint)cct__NextPowerOfTwo(baseFaceCount));
	}

	public static int cct__MinCbtDepth(CatmullClarkGPU.cc_SubdGPU subd)
	{
		int baseFaceCount = cct__BaseFaceCount(subd);

		return cct__FindMSB((uint)cct__NextPowerOfTwo(baseFaceCount));
	}


	/*******************************************************************************
	 * MaxBisectionDepth -- returns the maximum bisection depth supported by the subd
	 *
	 */
	public static int cct__MaxBisectionDepth(CatmullClark.cc_Subd subd)
	{
		return (CatmullClark.ccs_MaxDepth(subd) << 1) - 1;
	}

	public static int cct__MaxBisectionDepth(CatmullClarkGPU.cc_SubdGPU subd)
	{
		return (CatmullClarkGPU.ccs_MaxDepth(subd) << 1) - 1;
	}


	/*******************************************************************************
	 * GhostFaceCount -- returns the number of ghost nodes created in the CBT
	 *
	 */
	public static int cct__GhostFaceCount(CatmullClark.cc_Subd subd)
	{
		return (1 << cct__MinCbtDepth(subd)) - cct__BaseFaceCount(subd);
	}


	/*******************************************************************************
	 * BisectionDepth -- returns the bisection depth
	 *
	 */
	public static int cct_NodeBisectionDepth(ConcurrentBinaryTree.Node node, CatmullClark.cc_Subd subd)
	{
		int minDepth = cct__MinCbtDepth(subd);

		return node.depth - minDepth;
	}


	/*******************************************************************************
	 * Create -- Creates a CBT suitable for computing the tessellation (initializes at correct depth)
	 *
	 */
	public static ConcurrentBinaryTree cct_Create(CatmullClark.cc_Subd subd)
	{
		int minDepth = cct__MinCbtDepth(subd);
		int cbtDepth = minDepth + cct__MaxBisectionDepth(subd);

		return new ConcurrentBinaryTree(cbtDepth, minDepth);
	}

	public static ConcurrentBinaryTree cct_Create(CatmullClarkGPU.cc_SubdGPU subd)
	{
		int minDepth = cct__MinCbtDepth(subd);
		int cbtDepth = minDepth + cct__MaxBisectionDepth(subd);

		return new ConcurrentBinaryTree(cbtDepth, minDepth);
	}


	/*******************************************************************************
	 * FaceCount -- Returns the number of faces in the tessellation
	 *
	 */
	public static int cct_FaceCount(ConcurrentBinaryTree cbt, CatmullClark.cc_Subd subd)
	{
		int x = cbt.NodeCount();
		int y = cct__GhostFaceCount(subd);
		return cbt.NodeCount() - cct__GhostFaceCount(subd);
	}
}
