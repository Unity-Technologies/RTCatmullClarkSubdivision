/* public domain library for building binary trees in parallel
by Jonathan Dupuy
*/

#ifdef CCM_BF_FLAG_WRITE
RWStructuredBuffer<uint> u_Bitfield;
#else
StructuredBuffer<uint> u_Bitfield;
#endif


/*******************************************************************************
 * Size -- Returns the size of the bitfield
 * 
 */
int ccm_bf_Size()
{
	return u_Bitfield[0];
}


/*******************************************************************************
 * NextPowerOfTwo -- Returns the upper power of two value
 *
 * if the input is already a power of two, its value is returned.
 *
 */
int ccm_bf_NextPowerOfTwo(int x)
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
 * FindMSB -- Returns the position of the most significant bit
 *
 */
int ccm_bf_FindMSB(int x)
{
	int msb = 0;

	while (x > 1u)
	{
		++msb;
		x = x >> 1;
	}

	return msb;
}


/*******************************************************************************
 * SetBit -- Set a specific bit to either 0 or 1 in the bitfield
 *
 */
#ifdef CCM_BF_FLAG_WRITE
void ccm_bf_SetBit(int bitID, int bitValue)
{
	int offset = ccm_bf_Size();
	u_Bitfield[offset + bitID] = bitValue;
}
#endif


/*******************************************************************************
 * GetBit -- Get a specific bit in the bitfield
 *
 */
int ccm_bf_GetBit(int bitID)
{
	int offset = ccm_bf_Size();
	return u_Bitfield[offset + bitID];
}


/*******************************************************************************
 * BitCount -- Returns the number of bits set to one in the bit field
 *
 */
int ccm_bf_BitCount()
{
	return u_Bitfield[1];
}

/*******************************************************************************
 * DecodeNode -- Returns the leaf node associated to index nodeID
 *
 * This is procedure is for iterating over the one-valued bits.
 *
 */
int ccm_bf_DecodeBit(int handle)
{
	int bitID = 1;
	int bitFieldSize = ccm_bf_Size();

	while (bitID < bitFieldSize)
	{
		int heapValue = u_Bitfield[bitID * 2];
		int b = (int)handle < heapValue ? 0 : 1;

		bitID = 2 * bitID | b;
		handle -= heapValue * b;
	}

	return (bitID ^ bitFieldSize);
}

/*******************************************************************************
 * EncodeNode -- Returns the handle associated with the corresponding bitID
 *
 * This does the inverse of the DecodeNode routine. Note that this mapping
 * has the property that any bit set to 0 will be mapped to the ID of the next
 * bit set to one in the bit field.
 *
 */
int ccm_bf_EncodeBit(int bitID)
{
	int bitFieldSize = ccm_bf_Size();
	int arrayID = bitID + bitFieldSize;
	int handle = 0;

	while (arrayID > 1u)
	{
		uint sibling = (uint)arrayID & (~1u);
		int bitCount = u_Bitfield[sibling];

		handle += (int)(arrayID & 1u) * bitCount;
		arrayID = arrayID / 2;
	}

	return handle;
}