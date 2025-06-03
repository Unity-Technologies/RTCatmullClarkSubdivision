using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ConcurrentBinaryTree
{
	public uint[] m_buffer;


	/*******************************************************************************
	 * CBT Constructor
	 *
	 */
	public ConcurrentBinaryTree(int maxDepth, int initDepth = 0)
	{
		m_buffer = new uint[HeapUintSize(maxDepth)];
		m_buffer[0] = 1u << maxDepth;

		ResetToDepth(initDepth);
	}


	/*******************************************************************************
	 * ResetToDepth -- Resets a CBT to a specific subdivision level
	 *
	 */
	public void ResetToDepth(int depth)
	{
		uint minNodeID = 1u << depth;
		uint maxNodeID = 2u << depth;

		ClearBitField();

		for (uint nodeID = minNodeID; nodeID < maxNodeID; ++nodeID)
		{
			HeapWrite_BitField(new Node(nodeID, depth), 1u);
		}

		//ComputeSumReduction(); //XXX: disabled -- make sure to evaluate on the GPU !
	}
	public void ResetToRoot() { ResetToDepth(0); }
	public void ResetToCeil() { ResetToDepth(MaxDepth()); }


	/*******************************************************************************
	 * MaxDepth -- Returns the maximum depth of the CBT
	 *
	 */
	public int MaxDepth()
	{
		return FindLSB(m_buffer[0]);
	}

	/*******************************************************************************
	 * NodeCount -- Returns the number of nodes in the CBT
	 *
	 */
	public int NodeCount()
	{
		return (int)HeapRead(new Node(1u, 0));
	}


	/*******************************************************************************
	 * GetHeap -- Returns a read-only pointer to the heap memory
	 *
	 */
	public uint[] GetHeap()
	{
		return m_buffer;
	}


	/*******************************************************************************
	 * HeapByteSize -- Returns the amount of bytes consumed by the CBT heap
	 *
	 */
	public int HeapByteSize()
	{
		return HeapByteSize(MaxDepth());
	}


	/*******************************************************************************
	 * HeapUintSize -- Returns the number of uints consumed by the CBT heap
	 *
	 */
	public int HeapUintSize()
	{
		return HeapUintSize(MaxDepth());
	}


	/*******************************************************************************
	 * DecodeNode -- Returns the leaf node associated to index nodeID
	 *
	 * This is procedure is for iterating over the nodes.
	 *
	 */
	public Node DecodeNode(int handle)
	{
		Debug.Assert(handle < NodeCount(), "handle > NodeCount");
		Debug.Assert(handle >= 0, "handle < 0");

		Node node = new Node(1u, 0);

		while (HeapRead(node) > 1u)
		{
			Node heapNode = new Node(node.id <<= 1, ++node.depth);
			uint cmp = HeapRead(heapNode);
			uint b = (uint)handle < cmp ? 0u : 1u;

			node.id |= b;
			handle -= (int)(cmp * b);
		}

		return node;
	}


	/*******************************************************************************
	 * EncodeNode -- Returns the bit index associated with the Node
	 *
	 * This does the inverse of the DecodeNode routine.
	 *
	 */
	public int EncodeNode(Node node)
	{
		Debug.Assert(IsLeafNode(node), "node is not a leaf");

		int handle = 0;
		Node nodeIterator = node;

		while (nodeIterator.id > 1u)
		{
			Node sibling = LeftSiblingNode_Fast(nodeIterator);
			uint nodeCount = HeapRead(sibling);

			handle += (int)((nodeIterator.id & 1u) * nodeCount);
			nodeIterator = ParentNode_Fast(nodeIterator);
		}

		return handle;
	}


	/*******************************************************************************
	 * ComputeSumReduction -- Builds the sum reduction tree
	 *
	 */
	private void ComputeSumReduction()
	{
		for (int depth = MaxDepth() - 1; depth >= 0; --depth)
		{
			uint minNodeID = 1u << depth;
			uint maxNodeID = 2u << depth;

			for (uint nodeID = minNodeID; nodeID < maxNodeID; ++nodeID)
			{
				uint x0 = HeapRead(new Node(nodeID << 1, depth + 1));
				uint x1 = HeapRead(new Node(nodeID << 1 | 1, depth + 1));

				HeapWrite(new Node(nodeID, depth), x0 + x1);
			}
		}
	}


	/*******************************************************************************
	 * HeapByteSize -- Computes the number of Bytes to allocate for the CBT heap
	 *
	 * For a tree of max depth D, the number of Bytes is 2^(D-1).
	 *
	 */
	static public int HeapByteSize(int maxDepth)
	{
		return 1 << (maxDepth - 1);
	}


	/*******************************************************************************
	 * HeapUintSize -- Computes the number of uints to allocate for the bitfield
	 *
	 */
	static private int HeapUintSize(int maxDepth)
	{
		return HeapByteSize(maxDepth) >> 2;
	}


	/*******************************************************************************
	 * FindLSB -- Returns the position of the least significant bit
	 *
	 */
	static private int FindLSB(uint x)
	{
		int lsb = 0;

		while (((x >> lsb) & 1u) == 0u)
		{
			++lsb;
		}

		return lsb;
	}


	/*******************************************************************************
	 * MinValue -- Returns the minimum value between two inputs
	 *
	 */
	static private int MinValue(int a, int b)
	{
		return a < b ? a : b;
	}


	/*******************************************************************************
	 * ClearHeap -- Clears the heap
	 *
	 */
	private void ClearBitField()
	{
		int maxDepth = MaxDepth();
		int bufferMinID = 1 << (maxDepth - 4);
		int bufferMaxID = HeapUintSize(maxDepth);

		// TODO: FIX
		if (bufferMinID < 0)
			bufferMinID = 0;

		for (int bufferID = bufferMinID; bufferID < bufferMaxID; ++bufferID)
		{
			if(bufferID < 0 || bufferID >= m_buffer.Length)
			{
				int x = 0;
			}
			m_buffer[bufferID] = 0;
		}
	}


	/*******************************************************************************
	 * BitfieldInsert -- Inserts bitData in range [bitOffset, bitOffset + bitCount - 1]
	 *
	 */
	static private void
	BitFieldInsert(ref uint bitField, int bitOffset, int bitCount, uint bitData)
	{
		Debug.Assert(bitOffset < 32 && bitCount <= 32 && bitOffset + bitCount <= 32);
		uint bitMask = ~(~(0xFFFFFFFFu << bitCount) << bitOffset);

		bitField &= bitMask;
		bitField |= (bitData << bitOffset);
	}


	/*******************************************************************************
	 * BitFieldExtract -- Extracts bits [bitOffset, bitOffset + bitCount - 1] from
	 * a bitfield, returning them in the least significant bits of the result.
	 *
	 */
	static private uint
	BitFieldExtract(uint bitField, int bitOffset, int bitCount)
	{
		Debug.Assert(bitOffset < 32 && bitCount < 32 && bitOffset + bitCount <= 32);
		uint bitMask = ~(0xFFFFFFFFu << bitCount);

		return (bitField >> bitOffset) & bitMask;
	}


	/*******************************************************************************
	 * SetBitValue -- Sets the value of a bit stored in a bitfield
	 *
	 */
	static private void SetBitValue(ref uint bitField, int bitID, uint bitValue)
	{
		uint bitMask = ~(1u << bitID);

		bitField &= bitMask;
		bitField |= (bitValue << bitID);
	}


	/*******************************************************************************
	 * HeapWrite_BitField -- Sets the bit associated to a leaf node to bitValue
	 *
	 * This is a dedicated routine to write directly to the bitfield.
	 *
	 */
	private void HeapWrite_BitField(Node node, uint bitValue)
	{
		int bitID = NodeBitID_BitField(node);

		SetBitValue(ref m_buffer[bitID >> 5], bitID, bitValue);
	}


	/*******************************************************************************
	 * HeapArgs
	 *
	 * The LEB heap data structure uses an array of 32-bit words to store its data.
	 * Whenever we need to access a certain bit range, we need to query two such
	 * words (because sometimes the requested bit range overlaps two 32-bit words).
	 * The HeapArg data structure provides arguments for reading from and/or
	 * writing to the two 32-bit words that bound the queries range.
	 *
	 */
	private struct HeapArgs
	{
		public int bufferIndexLSB;
		public int bufferIndexMSB;
		public int bitOffsetLSB;
		public int bitCountLSB, bitCountMSB;
	}

	static private HeapArgs
	CreateHeapArgs(ConcurrentBinaryTree cbt, Node node, int bitCount)
	{
		int alignedBitOffset = cbt.NodeBitID(node);
		int maxBufferIndex = HeapUintSize(cbt.MaxDepth()) - 1;
		HeapArgs args;

		args.bufferIndexLSB = (alignedBitOffset >> 5);
		args.bufferIndexMSB = MinValue(args.bufferIndexLSB + 1, maxBufferIndex);
		args.bitOffsetLSB = alignedBitOffset & 31;
		args.bitCountLSB = MinValue(32 - args.bitOffsetLSB, bitCount);
		args.bitCountMSB = bitCount - args.bitCountLSB;

		return args;
	}


	/*******************************************************************************
	 * HeapWrite -- Sets bitCount bits located at nodeID to bitData
	 *
	 * Note that this procedure writes to at most two uint64 elements.
	 * Two elements are relevant whenever the specified interval overflows 64-bit
	 * words.
	 *
	 */
	private void HeapWriteExplicit(Node node, int bitCount, uint bitData)
	{
		HeapArgs args = CreateHeapArgs(this, node, bitCount);

		BitFieldInsert(ref m_buffer[args.bufferIndexLSB],
					   args.bitOffsetLSB,
					   args.bitCountLSB,
					   bitData);
		BitFieldInsert(ref m_buffer[args.bufferIndexMSB],
					   0,
					   args.bitCountMSB,
					   bitData >> args.bitCountLSB);
	}

	private void HeapWrite(Node node, uint bitData)
	{
		HeapWriteExplicit(node, NodeBitSize(node), bitData);
	}


	/*******************************************************************************
	 * HeapRead -- Returns bitCount bits located at nodeID
	 *
	 * Note that this procedure reads from two uint64 elements.
	 * This is because the data is not necessarily aligned with 64-bit
	 * words.
	 *
	 */
	private uint HeapReadExplicit(Node node, int bitCount)
	{
		HeapArgs args = CreateHeapArgs(this, node, bitCount);
		uint lsb = BitFieldExtract(m_buffer[args.bufferIndexLSB],
								   args.bitOffsetLSB,
								   args.bitCountLSB);
		uint msb = BitFieldExtract(m_buffer[args.bufferIndexMSB],
								   0,
								   args.bitCountMSB);

		return (lsb | (msb << args.bitCountLSB));
	}

	public uint HeapRead(Node node)
	{
		return HeapReadExplicit(node, NodeBitSize(node));
	}


	/*******************************************************************************
	 * NodeBitID -- Returns the bit index that stores data associated with a given node
	 *
	 * For a tree of max depth D and given an index in [0, 2^(D+1) - 1], this
	 * functions is used to emulate the behaviour of a lookup in an array, i.e.,
	 * uint[nodeID]. It provides the first bit in memory that stores
	 * information associated with the element of index nodeID.
	 *
	 * For data located at level d, the bit offset is 2^d x (3 - d + D)
	 * We then offset this quantity by the index by (nodeID - 2^d) x (D + 1 - d)
	 * Note that the null index (nodeID = 0) is also supported.
	 *
	 */
	public int NodeBitID(Node node)
	{
		int tmp1 = 2 << node.depth;
		int tmp2 = 1 + MaxDepth() - node.depth;

		return tmp1 + (int)node.id * tmp2;
	}


	/*******************************************************************************
	 * NodeBitID_BitField -- Computes the bitfield bit location associated to a node
	 *
	 * Here, the node is converted into a final node and its bit offset is
	 * returned, which is finalNodeID + 2^{D + 1}
	 */
	private int NodeBitID_BitField(Node node)
	{
		return NodeBitID(CeilNode(node));
	}


	/*******************************************************************************
	 * NodeBitSize -- Returns the number of bits storing the input node value
	 *
	 */
	private int NodeBitSize(Node node)
	{
		return MaxDepth() - node.depth + 1;
	}


	/*******************************************************************************
	 * Node Data Type
	 *
	 */
	public struct Node
	{
		public uint id;
		public int depth;

		public Node(uint nodeID, int nodeDepth)
		{
			id = nodeID;
			depth = nodeDepth;
		}
	}


	/*******************************************************************************
	 * IsCeilNode -- Checks if a node is a ceil node, i.e., that can not split further
	 *
	 */
	public bool IsCeilNode(Node node)
	{
		return (node.depth == MaxDepth());
	}


	/*******************************************************************************
	 * IsRootNode -- Checks if a node is a root node
	 *
	 */
	public bool IsRootNode(Node node)
	{
		return (node.id == 1u);
	}


	/*******************************************************************************
	 * IsNullNode -- Checks if a node is a null node
	 *
	 */
	public bool IsNullNode(Node node)
	{
		return (node.id == 0u);
	}


	/*******************************************************************************
	 * IsLeafNode -- Checks if a node is a leaf node, i.e., that has no descendants
	 *
	 */
	public bool IsLeafNode(Node node)
	{
		return (HeapRead(node) == 1u);
	}


	/*******************************************************************************
	 * ParentNode -- Computes the parent of the input node
	 *
	 */
	public Node ParentNode_Fast(Node node)
	{
		return new Node(node.id >> 1, node.depth - 1);
	}

	public Node cbt_ParentNode(Node node)
	{
		return IsNullNode(node) ? node : ParentNode_Fast(node);
	}


	/*******************************************************************************
	 * CeilNode -- Returns the associated ceil node, i.e., the deepest possible leaf
	 *
	 */
	public Node CeilNode_Fast(Node node)
	{
		int maxDepth = MaxDepth();

		return new Node(node.id << (maxDepth - node.depth), maxDepth);
	}

	public Node CeilNode(Node node)
	{
		return IsNullNode(node) ? node : CeilNode_Fast(node);
	}


	/*******************************************************************************
	 * SiblingNode -- Computes the sibling of the input node
	 *
	 */
	public Node SiblingNode_Fast(Node node)
	{
		return new Node(node.id ^ 1u, node.depth);
	}

	public Node SiblingNode(Node node)
	{
		return IsNullNode(node) ? node : SiblingNode_Fast(node);
	}


	/*******************************************************************************
	 * RightSiblingNode -- Computes the right sibling of the input node
	 *
	 */
	public Node RightSiblingNode_Fast(Node node)
	{
		return new Node(node.id | 1u, node.depth);
	}

	public Node RightSiblingNode(Node node)
	{
		return IsNullNode(node) ? node : RightSiblingNode_Fast(node);
	}


	/*******************************************************************************
	 * LeftSiblingNode -- Computes the left sibling of the input node
	 *
	 */
	public Node LeftSiblingNode_Fast(Node node)
	{
		return new Node(node.id & (~1u), node.depth);
	}

	public Node LeftSiblingNode(Node node)
	{
		return IsNullNode(node) ? node : LeftSiblingNode_Fast(node);
	}


	/*******************************************************************************
	 * RightChildNode -- Computes the right child of the input node
	 *
	 */
	public Node RightChildNode_Fast(Node node)
	{
		return new Node((node.id << 1) | 1u, node.depth + 1);
	}

	public Node RightChildNode(Node node)
	{
		return IsNullNode(node) ? node : RightChildNode_Fast(node);
	}


	/*******************************************************************************
	 * LeftChildNode -- Computes the left child of the input node
	 *
	 */
	public Node LeftChildNode_Fast(Node node)
	{
		return new Node((node.id << 1), node.depth + 1);
	}

	public Node cbt_LeftChildNode(Node node)
	{
		return IsNullNode(node) ? node : LeftChildNode_Fast(node);
	}

}
