// accessors
int cct_FaceCount(const int cbtID, const int subdID);
int cct_NodeBisectionDepth(const in cbt_Node node, const int subdID);

// face attributes
struct cct_Face
{
	int3 halfedgeIDs;
};
cct_Face cct_DecodeFace(const in cbt_Node node, const int subdID);

// diamond parent
struct cct_DiamondParent
{
	cbt_Node base, top;
};
cct_DiamondParent cct_DecodeDiamondParent(const in cbt_Node node,
										  const int subdID);

// split / merge
#ifdef CBT_FLAG_WRITE
void cct_SplitNode(const int cbtID,
				   const in cbt_Node node,
				   const int subdID);
void cct_MergeNode(const int cbtID,
				   const in cbt_Node node,
				   const in cct_DiamondParent diamond,
				   const int subdID);
#endif


/*******************************************************************************
 * NextPowerOfTwo -- Returns the closest power of two
 *
 */
int cct__NextPowerOfTwo(int x)
{
	x--;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	x++;

	return x;
}


/*******************************************************************************
 * BaseFaceCount -- returns the number of faces at root
 *
 */
int cct__BaseFaceCount(const int subdID)
{
	return ccm_HalfedgeCount;
}


/*******************************************************************************
 * MinCbtDepth -- returns the minimum depth required for the CBT to work
 *
 */
int cct__MinCbtDepth(const int subdID)
{
	const int baseFaceCount = cct__BaseFaceCount(subdID);

	return firstbithigh(cct__NextPowerOfTwo(baseFaceCount));
}


/*******************************************************************************
 * MaxBisectionDepth -- returns the maximum bisection depth supported by the subd
 *
 */
int cct__MaxBisectionDepth(const int subdID)
{
	return ((ccs_MaxDepth(subdID) << 1) - 1);
}


/*******************************************************************************
 * GhostFaceCount -- returns the number of ghost nodes created in the CBT
 *
 */
int cct__GhostFaceCount(const int subdID)
{
	return (1 << cct__MinCbtDepth(subdID)) - cct__BaseFaceCount(subdID);
}


/*******************************************************************************
 * BisectionDepth -- returns the bisection depth
 *
 */
int cct_NodeBisectionDepth(const in cbt_Node node, const int subdID)
{
	const int minDepth = cct__MinCbtDepth(subdID);

	return node.depth - minDepth;
}


/*******************************************************************************
 * FaceCount -- Returns the number of faces in the tessellation
 *
 */
int cct_FaceCount(const int cbtID, const int subdID)
{
	return int(cbt_NodeCount() - cct__GhostFaceCount(subdID));
}


/*******************************************************************************
 * NeighborFaceIDs
 *
 */
struct cct__NeighborFaceIDs
{
	int left, right, edge, face;
};

cct__NeighborFaceIDs cct__CreateNeighborFaceIDs(int left, int right, int edge, int face)
{
	cct__NeighborFaceIDs faceIDs;

	faceIDs.left = left;
	faceIDs.right = right;
	faceIDs.edge = edge;
	faceIDs.face = face;

	return faceIDs;
}


/*******************************************************************************
 * CreateDiamondParent -- Private DiamondParent constructor
 *
 */
cct_DiamondParent cct__CreateDiamondParent(const in cbt_Node base, const in cbt_Node top)
{
	cct_DiamondParent diamond;

	diamond.base = base;
	diamond.top = top;

	return diamond;
}


/*******************************************************************************
 * GetBitValue -- Returns the value of a bit stored in a 32-bit word
 *
 */
uint cct__GetBitValue(const uint bitField, int bitID)
{
	return ((bitField >> bitID) & 1u);
}


/*******************************************************************************
 * EvenRule -- Even half-edge splitting rule
 *
 */
void cct__EvenRule(inout int3 x, uint b)
{
	if (b == 0)
	{
		x[2] = (x[0] | 2);
		x[1] = (x[0] | 1);
		x[0] = (x[0] | 0);
	}
	else
	{
		x[0] = (x[2] | 2);
		x[1] = (x[2] | 3);
		x[2] = (x[2] | 0);
	}
}


/*******************************************************************************
 * EvenRule -- Odd half-edge splitting rule
 *
 */
void cct__OddRule(inout int3 x, uint b)
{
	if (b == 0)
	{
		x[2] = (x[1] << 2) | 0;
		x[1] = (x[0] << 2) | 2;
		x[0] = (x[0] << 2) | 0;
	}
	else
	{
		x[0] = (x[1] << 2) | 0;
		x[1] = (x[1] << 2) | 2;
		x[2] = (x[2] << 2) | 0;
	}
}


/*******************************************************************************
 * CreateFace -- Face constructor
 *
 */
cct_Face cct__CreateFace(int h0, int h1, int h2)
{
	cct_Face face;

	face.halfedgeIDs[0] = h0;
	face.halfedgeIDs[1] = h1;
	face.halfedgeIDs[2] = h2;

	return face;
}


/*******************************************************************************
 * BaseHalfedgeID -- Decodes
 *
 */
int cct__BaseHalfedgeID(const cbt_Node node, const int subdID)
{
	const int bisectionDepth = cct_NodeBisectionDepth(node, subdID);
	const int faceID = int(node.id ^ (1 << node.depth));
	const int baseHalfedgeID = faceID >> bisectionDepth;

	return baseHalfedgeID;
}


/*******************************************************************************
 * DecodeBaseFace -- Decodes
 *
 */
cct_Face cct__BaseFace(const in cbt_Node node, const int subdID)
{
	const int baseHalfedgeID = cct__BaseHalfedgeID(node, subdID);

	return cct__CreateFace
	(
		4 * baseHalfedgeID,
		4 * baseHalfedgeID + 2,
		4 * ccm_HalfedgeNextID(subdID, baseHalfedgeID)
	);
}


/*******************************************************************************
 * DecodeFace -- Face constructor
 *
 */
cct_Face cct_DecodeFace(const in cbt_Node node, const int subdID)
{
	int bisectionDepth = cct_NodeBisectionDepth(node, subdID);
	cct_Face face = cct__BaseFace(node, subdID);
	int pingPong = 0;

	for (int bitID = bisectionDepth - 1; bitID >= 0; --bitID)
	{
		const uint bitValue = cct__GetBitValue(node.id, bitID);

		if ((pingPong & 1) == 0)
		{
			cct__EvenRule(face.halfedgeIDs, bitValue);
		}
		else
		{
			cct__OddRule(face.halfedgeIDs, bitValue);
		}
		pingPong ^= 1;
	}

	// swap winding
	if ((bisectionDepth & 1) == 0)
	{
		const int tmp = face.halfedgeIDs[0];

		face.halfedgeIDs[0] = face.halfedgeIDs[2];
		face.halfedgeIDs[2] = tmp;
	}

	return face;
}


/*******************************************************************************
 * LoadNeighborFaceIDs -- Loads root neighborhood info
 *
 */
cct__NeighborFaceIDs cct__BaseNeighborFaceIDs(const in cbt_Node node, const int subdID)
{
#if 0
    const cct_Face face = cct__BaseFace(node, subdID);
    const int minDepth = cct__MinCbtDepth(subdID);
    const int n0 = ccs_HalfedgeTwinID(subdID, face.halfedgeIDs[0], 1);
    const int n1 = ccs_HalfedgeTwinID(subdID, face.halfedgeIDs[1], 1);
    const int leftID  = (1 << minDepth) | (n0 >> 1);
    const int rightID = (1 << minDepth) | (n1 >> 1);
    const int faceID = int(node.id >> cct_NodeBisectionDepth(node, subdID));
    const int edgeID = faceID ^ 1;
#else
	const int baseHalfedgeID = cct__BaseHalfedgeID(node, subdID);
	const int twinID = ccm_HalfedgeTwinID(subdID, baseHalfedgeID);
	const int prevID = ccm_HalfedgeNextID(subdID, baseHalfedgeID);
	const int nextID = ccm_HalfedgePrevID(subdID, baseHalfedgeID);
	const int minDepth = cct__MinCbtDepth(subdID);
	const int bitMask = 1 << minDepth;
	const int leftID = bitMask | nextID;
	const int rightID = bitMask | prevID;
	const int edgeID = bitMask | twinID;
	const int faceID = bitMask | baseHalfedgeID;

#endif

	return cct__CreateNeighborFaceIDs(leftID, rightID, edgeID, faceID);
}


/*******************************************************************************
 * SubdivideNeighborFaceIDs -- Computes new neighborhood after one subdivision step
 *
 */
cct__NeighborFaceIDs cct__SubdivideNeighborFaceIDs(const in cct__NeighborFaceIDs neighborFaceIDs, uint bitValue)
{
	const int L = neighborFaceIDs.left, R = neighborFaceIDs.right,
              E = neighborFaceIDs.edge, F = neighborFaceIDs.face;

	if (bitValue == 0)
	{
		return cct__CreateNeighborFaceIDs(E << 1 | 1, F << 1 | 1, L << 1 | 1, F << 1 | 0);
	}
	else
	{
		return cct__CreateNeighborFaceIDs(F << 1 | 0, E << 1 | 0, R << 1 | 0, F << 1 | 1);
	}
}


/*******************************************************************************
 * SubdivideNeighborFaceIDs -- Computes new neighborhood after one subdivision step
 *
 */
cct__NeighborFaceIDs cct__DecodeNeighborFaceIDs(const in cbt_Node node, const int subdID)
{
	const int bisectionDepth = cct_NodeBisectionDepth(node, subdID);
	cct__NeighborFaceIDs neighborFaceIDs = cct__BaseNeighborFaceIDs(node, subdID);

	for (int bitID = bisectionDepth - 1; bitID >= 0; --bitID)
	{
		const uint bitValue = cct__GetBitValue(node.id, bitID);

		neighborFaceIDs = cct__SubdivideNeighborFaceIDs(neighborFaceIDs,
                                                        bitValue);
	}

	return neighborFaceIDs;
}


/*******************************************************************************
 * SplitNode -- Bisects a triangle in the current tessellation
 *
 */
#ifdef CBT_FLAG_WRITE
void cct_SplitNode(const int cbtID, const cbt_Node node, const int subdID)
{
	cbt_Node nodeIterator = node;

	cbt_SplitNode_Fast(nodeIterator);

	int safeGuard = 0;
	while (nodeIterator.id >= (1 << cct__MinCbtDepth(subdID)) && safeGuard < 100)
	{
		cct__NeighborFaceIDs faceIDs = cct__DecodeNeighborFaceIDs(nodeIterator, subdID);

		if (faceIDs.edge <= 0)
			break;

		nodeIterator.id = faceIDs.edge;
		cbt_SplitNode_Fast(nodeIterator);
		nodeIterator = cbt_ParentNode_Fast(nodeIterator);
		cbt_SplitNode_Fast(nodeIterator);

		safeGuard++;
	}
}
#endif


/*******************************************************************************
 * DecodeDiamondParent -- Decodes the diamond associated to the Node
 *
 * If the neighbour part does not exist, the parentNode is copied instead.
 *
 */
cct_DiamondParent cct_DecodeDiamondParent(in const cbt_Node node, const int subdID)
{
	cbt_Node parentNode = cbt_ParentNode_Fast(node);
	int edgeTwinID = cct__DecodeNeighborFaceIDs(parentNode, subdID).edge;
	cbt_Node edgeNeighborNode = cbt_CreateNode
	(
		edgeTwinID > 0 ? edgeTwinID : parentNode.id,
		parentNode.depth
	);

	return cct__CreateDiamondParent(parentNode, edgeNeighborNode);
}


/*******************************************************************************
 * HasDiamondParent -- Determines whether a diamond parent is actually stored
 *
 * This procedure checks that the diamond parent is encoded in the CBT.
 * We can perform this test by checking that both the base and top nodes
 * that form the diamond parent are split, i.e., CBT[base] = CBT[top] = 2.
 * This is a crucial operation for implementing the leb_Merge routine.
 *
 */
bool cct__HasDiamondParent(const int cbtID, in const cct_DiamondParent diamondParent)
{
	bool canMergeBase = cbt_HeapRead(diamondParent.base) <= 2u;
	bool canMergeTop = cbt_HeapRead(diamondParent.top) <= 2u;

	return canMergeBase && canMergeTop;
}


/*******************************************************************************
 * MergeNode -- Bisects a triangle in the current tessellation
 *
 */
bool cct_IsRootNode(in const cbt_Node node, const int subdID)
{
	return node.depth == cct__MinCbtDepth(subdID);
}

#ifdef CBT_FLAG_WRITE
void cct_MergeNode(const int cbtID, in const cbt_Node node, in const cct_DiamondParent diamondParent, const int subdID)
{
#if 0
	if (!cct_IsRootNode(node, subdID) && cct__HasDiamondParent(cbtID, diamondParent))
	{
		cbt_MergeNode(node);
	}
#else
    const cbt_Node sibling = cbt_SiblingNode_Fast(node);
    const cbt_Node left = cbt_LeftChildNode_Fast(diamondParent.top);
    const cbt_Node right = cbt_RightChildNode_Fast(diamondParent.top);
    bool leafTest = cbt_IsLeafNode(sibling)
                  && cbt_IsLeafNode(left)
                  && cbt_IsLeafNode(right);

    if (!cct_IsRootNode(node, subdID) && leafTest) {
        cbt_MergeNode(node);
        // FIXME
        cbt_MergeNode_Fast(cbt_RightChildNode_Fast(diamondParent.top));
    }
#endif
}
#endif
