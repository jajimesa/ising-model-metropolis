! Simulación de Monte Carlo del modelo de Ising en red cuadrada 2D
! con condiciones de contorno periódicas, constante de intercambio J=1
program Ising_2D_Metropolis

    implicit double precision(a-h,o-z)
    parameter (L=120, N=L*L, NT=60)  ! L: longitud de la red, N: número total de espines, NT: número de pasos de temperatura

    integer s(N), n1(N), n2(N), n3(N), n4(N)  ! Espines y vecinos
    dimension h(-4:4)                         ! Factores de Boltzmann-Gibbs

    integer mc, mc_non_critical_inicial, mc_critical, mc_increment_steps
    integer M, M_non_critical, M_critical

    data M0 /1000/                            ! M0: pasos de termalización

    data mc_non_critical_inicial, mc_critical /10, 100/ ! Valores inicial y crítico de mc
    mc_increment_steps = 5 ! Número de pasos para aumentar/disminuir mc

    data M_non_critical, M_critical /4096, 8192/ ! Número de muestras en regiones no críticas y crítica

    open(unit=77, file='fort.77', status='replace', action='write')
    print *, 'Archivo fort.77 abierto correctamente'

    ! Generar los vectores de los vecinos con condiciones de contorno periódicas
    call neighbors(n1, n2, n3, n4, L)

    ! Inicialización aleatoria de la configuración de espines
    do i = 1, N
        if (ran_u() < 0.5d0) then
            s(i) = +1
        else
            s(i) = -1
        endif
    enddo

    ! Barrido adaptativo de temperaturas en el intervalo [2.2, 2.36]
    do 999 iT = 0, NT - 1

        print *, 'It. ', iT

        ! Determinar la temperatura de acuerdo a los intervalos definidos
        if (iT .lt. 5) then
            ! Primer intervalo: [2.2, 2.25] con 5 pasos
            T = 2.2d0 + (2.25d0 - 2.2d0) * dble(iT) / 4.0d0
            mc = mc_non_critical_inicial + int(dble(mc_critical - mc_non_critical_inicial) * dble(iT) / dble(mc_increment_steps))
            M = M_non_critical + int(dble(M_critical - M_non_critical) * dble(iT) / dble(mc_increment_steps))
        elseif (iT .lt. NT - 5) then
            ! Intervalo crítico: [2.25, 2.32] con NT - 10 pasos
            T = 2.25d0 + (2.32d0 - 2.25d0) * dble(iT - 5) / dble(NT - 11)
            mc = mc_critical
            M = M_critical
        else
            ! Último intervalo: [2.32, 2.36] con 5 pasos
            T = 2.32d0 + (2.36d0 - 2.32d0) * dble(iT - (NT - 5)) / 4.0d0
            mc_reduction_step = (NT - 1 - iT)
            mc = mc_non_critical_inicial + int( &
                 dble(mc_critical - mc_non_critical_inicial) * &
                 dble(mc_reduction_step) / &
                 dble(mc_increment_steps) &
                 )
            M = M_non_critical + int( &
                dble(M_critical - M_non_critical) * &
                dble(mc_reduction_step) / &
                dble(mc_increment_steps) &
                )
        endif
        ! Asegurar que mc no sea menor que el valor inicial
        mc = max(mc, mc_non_critical_inicial)

        ! Crear el array con los factores de Boltzmann-Gibbs para cada temperatura
        do j = -4, 4, 2
            h(j) = min(1.0d0, exp(-2 * j / T))
        enddo

        ! Termalización: Relajación al estado estacionario
        do ij = 1, M0 * N

            ! if (mod(ij, 10000) == 0) print *, 'Termalización en paso ', ij

            i = i_ran(N)  ! Elegir un espín aleatorio
            ib = s(i) * (s(n1(i)) + s(n2(i)) + s(n3(i)) + s(n4(i)))  ! Energía local
            if (ran_u() < h(ib)) s(i) = -s(i)  ! Actualización según Metropolis
        enddo

        ! Inicialización de promedios
        c = 0.0d0
        rm = 0.0d0
        rm2 = 0.0d0

        ! Magnetización inicial (absoluta)
        rm1 = real(abs(sum(s))) / N

        ! Actualización y medición
        do im = 1, M
            do ij = 1, mc * N ! Usamos el número de pasos MC adaptado

                ! if (mod(ij, 100000) == 0) print *, 'Medición en im = ', im, ', ij = ', ij

                i = i_ran(N)  ! Seleccionar un espín aleatorio
                ib = s(i) * (s(n1(i)) + s(n2(i)) + s(n3(i)) + s(n4(i)))
                if (ran_u() < h(ib)) s(i) = -s(i)  ! Actualizar el espín
            enddo

            ! Medir magnetización
            rm0 = real(abs(sum(s))) / N

            ! Acumular valores para los promedios
            rm = rm + rm0
            rm2 = rm2 + rm0 * rm0
            c = c + rm0 * rm1
            rm1 = rm0
        enddo

        ! Calcular promedios finales
        rm = rm / M
        rm2 = rm2 / M - rm * rm
        c = (c / M - rm * rm) / rm2

        ! Calcular tiempo de autocorrelación y error
        if (c .ne. 1.0d0) then
            tau = c / (1.0d0 - c)
        else
            tau = 0.0d0
        endif
        error = sqrt(rm2 * (2 * tau + 1) / M)

        ! Escribir resultados: T, magnetización media, susceptibilidad, error, tiempo de autocorrelación, coeficiente de correlación
        write(77, '(f10.6, 1p5e16.6)') T, rm, rm2, error, mc * tau, c

        print *, 'Guardados datos de la iteración ', iT

999 continue

close(77)

end program Ising_2D_Metropolis


! Subrutina para calcular y guardar los vecinos más cercanos de un modelo de Ising 2D
! con condiciones de contorno periódicas
subroutine neighbors(n1,n2,n3,n4,L)
    dimension n1(L*L),n2(L*L),n3(L*L),n4(L*L)
    do ix=1,L
    do iy=1,L
    i=(iy-1)*L+ix
    ix1=ix+1
    if (ix1.eq.L+1) ix1=1
    n1(i)=(iy-1)*L+ix1
    iy2=iy+1
    if (iy2.eq.L+1) iy2=1
    n2(i)=(iy2-1)*L+ix
    ix3=ix-1
    if (ix3.eq.0) ix3=L
    n3(i)=(iy-1)*L+ix3
    iy4=iy-1
    if (iy4.eq.0) iy4=L
    n4(i)=(iy4-1)*L+ix
    enddo
    enddo
end subroutine neighbors


! Función para generar un índice aleatorio entre 1 y N
function i_ran(N)
    implicit none
    integer, intent(in) :: N
    integer :: i_ran
    real :: r

    ! Generar un número aleatorio en el rango [0, 1)
    call random_number(r)

    ! Escalarlo para que esté entre 1 y N
    i_ran = 1 + int(r * real(N))
end function i_ran


! Función para generar un número aleatorio uniforme entre 0 y 1
function ran_u() result(r)
    implicit none
    real(8) :: r  ! Cambiar el tipo a REAL(8)

    ! Generar un número aleatorio en el rango [0, 1)
    call random_number(r)
end function ran_u