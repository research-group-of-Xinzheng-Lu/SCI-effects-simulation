subroutine GenFile

use GenModelPara
USE IFPORT

implicit none
integer*4 ::tmp
inquire(file="FILES_MPI",exist=exist_f)
if(.not.exist_f) then
    tmp=MAKEDIRQQ("./FILES_MPI")
end if
inquire(file="MONITOR",exist=exist_f)
if(.not.exist_f) then
    tmp=MAKEDIRQQ("./MONITOR")
end if

return

end subroutine GenFile