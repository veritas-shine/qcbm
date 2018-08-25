import {ops, meta} from 'projectq'
const {Rx, Ry, Rz, CNOT} = ops
/**
 * @class CircuitBlock
 * the building block of a circuit. This is an abstract class.
 */
class CircuitBlock {
  /**
   * @constructor
   * @param {number} numBit
   */
  constructor(numBit) {
    this.numBit = numBit
  }

  /**
   * build a quantum circuit.

   Args:
   theta_list (1darray<float>, len=3*num*bit*(depth+1)): parameters in this quantum circuit, here, depth equals to the number of entanglement operations.

   Return:
   remaining theta_list
   * @param qureg
   * @param theta_list
   */
  call(qureg, theta_list) {

  }

  /**
   * number of parameters it consume.
   */
  get numParam() {

  }

  /**
   * build this block into a sequence of csr_matrices.

   Args:
   theta_list (1darray): parameters,

   Returns:
   list: a list of csr_matrices, apply them on a vector to perform operation.
   * @param theta_list
   */
  tocsr(theta_list) {

  }
}

/**
 * BlockQueue is a sequence of CircuitBlock instances.
 * @class BlockQueue
 */
class BlockQueue extends Array {

  get numBit() {
    return this[0].numBit
  }

  get numParam() {
    return this.reduce((pre, current) => pre + current.numParam, 0)
  }

  call(qureg, theta_list) {
    this.forEach(block => {
      const thetaI = theta_list.slice(0, block.numParam)
      theta_list = theta_list.slice(block.numParam)
      block.call(qureg, thetaI)
    })

    assert(theta_list.length === 0)
  }

  toString() {
    this.map(b => b.toString()).join('\n')
  }
}

/**
 * Clever Block Queue that keep track of theta_list changing history, for fast update.
 * @class CleverBlockQueue
 */
class CleverBlockQueue  extends BlockQueue {
  constructor() {
    super()
    this.theta_last = undefined
    this.memo = undefined
  }

  call(qureg, theta_list) {
    if ! isinstance(qureg, np.ndarray):
    return super(CleverBlockQueue, self).__call__(qureg, theta_list)
// cache? if theta_list change <= 1 parameters, then don't touch memory.
    remember = this.theta_last is None or (abs(this.theta_last-theta_list)>1e-12).sum() > 1

    mats = []
    theta_last = this.theta_last
    if remember:
    this.theta_last = theta_list.copy()

    qureg_ = qureg
    for iblock, block in enumerate(self):
// generate or use a block matrix
    numParam = block.numParam
    theta_i, theta_list = np.split(theta_list, [numParam])
    if theta_last is ! undefined:
      theta_o, theta_last = np.split(theta_last, [numParam])
    if this.memo is ! undefined && (numParam==0 or np.abs(theta_i-theta_o).max()<1e-12):
// use data cached in memory
    mat = this.memo[iblock]
  else:
    if this.memo is ! undefined && ! remember:
      // update the changed gate, but not touching memory.
      mat = _rot_tocsr_update1(block, this.memo[iblock], theta_o, theta_i)
  else:
// regenerate one
    mat = block.tocsr(theta_i)
    for mat_i in mat:
    qureg_ = mat_i.dot(qureg_)
    mats.append(mat)

    if remember:
// cache data
    this.memo = mats
// update register
    qureg[...] = qureg_
    assert(len(theta_list)==0)
  }
}

class ArbituaryRotation extends CircuitBlock {
  constructor(numBit) {
    super(numBit)
    /**
     * @type {boolean[]}
     */
    this.mask = new Array(3 * numBit)
    this.mask.fill(true)
  }

  call(qureg, theta_list) {
    const gates = [Rz, Rx, Rz]
    const theta_list_ = new Array(this.numBit * 3)
    theta_list_.fill(0)
    theta_list_[this.mask] = theta_list

    theta_list_.forEach((theta, i) => {
      const mask = this.mask[i]
      const ibit = Math.floor(i / 3)
      const igate = i % 3
      if (mask) {
        const gate = new gates[igate](theta)
        gate.or(qureg[ibit])
      }
    })
  }

  toString() {
    return `Rotate[${this.numParam}]`
  }

  get numParam() {
    return this.mask.reduce((pre, current) => pre + current, 0)
  }

  /**
   * transform this block to csr_matrix.
   * @param theta_list
   */
  tocsr(theta_list) {
    const theta_list_ = new Array(3 * this.numBit)
    theta_list_.fill(0)
    theta_list_[this.mask] = theta_list

    rots = [qclibs.rot(*ths) for ths in theta_list_.reshape([this.numBit,3])]
    res = [qclibs._([rot], [i], this.numBit) for i,rot in enumerate(rots)]
    return res
  }
}

class CNOTEntangler extends CircuitBlock {
  constructor(numBit, pairs) {
    super(numBit)
    this.pairs = pairs
  }

  toString() {
    const str = this.pairs.map((j, i) => `${j}-${i}`).join(',')
    return `CNOT(${str})`
  }

  call(qureg) {
    this.pairs.forEach(pair => CNOT.or(tuple(qureg[pair[0]], qureg[pair[1]])))
  }

  get numParam() {
    return 0
  }

  /**
   * transform this block to csr_matrix.
   * @param theta_list
   */
  tocsr(theta_list) {
    i, j = this.pairs[0]
    res = qclibs.CNOT(i, j, this.numBit)
    for i, j in this.pairs[1:]:
    res = qclibs.CNOT(i,j,this.numBit).dot(res)
    res.eliminate_zeros()
    return [res]
  }

  /**
   * rotation layer csr_matrix update method.

   Args:
   rot (ArbituaryRotation): rotatio layer.
   old (csr_matrix): old matrices.
   theta_old (1darray): old parameters.
   theta_new (1darray): new parameters.

   Returns:
   csr_matrix: new rotation matrices after the theta changed.
   * @param rot
   * @param old
   * @param theta_old
   * @param theta_new
   * @private
   */
  _rot_tocsr_update1(rot, old, theta_old, theta_new) {
    idiff_param = np.where(abs(theta_old-theta_new)>1e-12)[0].item()
    idiff = np.where(rot.mask)[0][idiff_param]

// get rotation parameters
    isite = idiff//3
    theta_list_ = np.zeros(3*rot.numBit)
    theta_list_[rot.mask] = theta_new

    new = old[:]
    new[isite] = qclibs._(qclibs.rot(*theta_list_[isite*3:isite*3+3]), isite, rot.numBit)
    return new
  }

  get_demo_circuit(numBit, depth, pairs) {
    blocks = []
// build circuit
    for idepth in range(depth+1):
    blocks.append(ArbituaryRotation(numBit))
    if idepth!=depth:
    blocks.append(CNOTEntangler(numBit, pairs))

// set leading and trailing Rz to disabled
    blocks[0].mask[::3] = False
    blocks[-1].mask[2::3] = False
    return CleverBlockQueue(blocks)
  }
}
