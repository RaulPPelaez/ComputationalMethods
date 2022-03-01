#include"vector.h"
namespace utils{
  inline real3 apply_mic(real3 pos, real3 boxSize){
    return pos - floorf((pos/boxSize + real(0.5)))*boxSize;
  }
  
  inline int3 indexToCell3D(int icell, int3 ncells){
    int3 cell;
    cell.x = icell%3 - 1;
    cell.y = (icell/3)%3 - 1;
    cell.z = (icell/9) - 1;
    return cell;
  }

  inline int3 pbc_cell(int3 cell, int3 ncells){
    if(cell.x < 0) cell.x += ncells.x;
    if(cell.x >= ncells.x) cell.x -= ncells.x;
    if(cell.y < 0) cell.y += ncells.y;
    if(cell.y >= ncells.y) cell.y -= ncells.y;
    if(cell.z < 0) cell.z += ncells.z;
    if(cell.z >= ncells.z) cell.z -= ncells.z;
    return cell;
  }

}
