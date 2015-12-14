module mod_matrix
    use mod_global
    implicit none

    interface tensor_prod
        module procedure tensor_prod_2,tensor_prod_3,   &
                         tensor_prod_4,tensor_prod_5,   &
                         tensor_prod_6
     end interface

    contains

    function tensor_prod_2(vec1,vec2)
    implicit none
    real(c_double),intent(in) :: vec1(DIM),vec2(DIM)
    real(c_double) :: tensor_prod_2(DIM,DIM)

    integer(c_int) :: p,q
    do q=1,DIM
        do p=1,DIM
            tensor_prod_2(p,q) = vec1(p)*vec2(q)
        enddo
    enddo
    return
    end function tensor_prod_2

    function tensor_prod_3(vec1,vec2,vec3)
    implicit none
    real(c_double),intent(in) :: vec1(DIM),vec2(DIM),vec3(DIM)
    real(c_double) :: tensor_prod_3(DIM,DIM,DIM)

    integer(c_int) :: p,q,r
    do r=1,DIM
        do q=1,DIM
            do p=1,DIM
                tensor_prod_3(p,q,r) = vec1(p)*vec2(q)*vec3(r)
            enddo
        enddo
    enddo
    return
    end function tensor_prod_3

    function tensor_prod_4(vec1,vec2,vec3,vec4)
    implicit none

    ! Passed variables
    real(c_double),intent(in) :: vec1(DIM),vec2(DIM),vec3(DIM),vec4(DIM)
    real(c_double)            :: tensor_prod_4(DIM,DIM,DIM,DIM)

    ! Local variables
    integer(c_int) :: p,q,r,s
    do s=1,DIM
        do r=1,DIM
            do q=1,DIM
                do p=1,DIM
                    tensor_prod_4(p,q,r,s) = vec1(p)*vec2(q)*vec3(r)*vec4(s)
                enddo
            enddo
        enddo
    enddo
    return
    end function tensor_prod_4

    function tensor_prod_5(vec1,vec2,vec3,vec4,vec5)
    implicit none
    real(c_double),intent(in) :: vec1(DIM),vec2(DIM),vec3(DIM),vec4(DIM),vec5(DIM)
    real(c_double) :: tensor_prod_5(DIM,DIM,DIM,DIM,DIM)

    integer(c_int) :: p,q,r,s,t
    do t=1,DIM
        do s=1,DIM
            do r=1,DIM
                do q=1,DIM
                    do p=1,DIM
                        tensor_prod_5(p,q,r,s,t) = vec1(p)*vec2(q)*vec3(r)*vec4(s)*vec5(t)
                    enddo
                enddo
            enddo
        enddo
    enddo
    return
    end function tensor_prod_5

    function tensor_prod_6(vec1,vec2,vec3,vec4,vec5,vec6)
    implicit none
    real(c_double),intent(in) :: vec1(DIM),vec2(DIM),vec3(DIM),vec4(DIM),vec5(DIM),vec6(DIM)
    real(c_double) :: tensor_prod_6(DIM,DIM,DIM,DIM,DIM,DIM)

    integer(c_int) :: p,q,r,s,t,u
    do u=1,DIM
        do t=1,DIM
            do s=1,DIM
                do r=1,DIM
                    do q=1,DIM
                        do p=1,DIM
                            tensor_prod_6(p,q,r,s,t,u) = vec1(p)*vec2(q)*vec3(r)*vec4(s)*vec5(t)*vec6(u)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    return
    end function tensor_prod_6

end module mod_matrix
