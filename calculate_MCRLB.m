function MCRLB = calculate_MCRLB(pseudotrue_estimated_states,f,f1,f2,f3,f4,h,h1,h2,Q,R,Rbar,true_measurements_in_polar,simulation_length,state_dimension,x0,P0,df_dx,dh_dx,df2_dx2,d2h1_dx2,d2h2_dx2,apply_numerical_derivative,apply_tridiagonal_method)                                 

Amatrix = zeros(state_dimension*simulation_length);
Bmatrix = zeros(state_dimension*simulation_length);
bvector = zeros(state_dimension*simulation_length,1);
invR=inv(R);
invQ=inv(Q);
for k = 1:simulation_length
    if apply_numerical_derivative
        % calculate the derivatives numerically
        first_derivative_of_f = d_dx(f,pseudotrue_estimated_states(:,k));

        second_derivative_of_f1 = d2_dx2(f1,pseudotrue_estimated_states(:,k));
        second_derivative_of_f2 = d2_dx2(f2,pseudotrue_estimated_states(:,k));
        second_derivative_of_f3 = d2_dx2(f3,pseudotrue_estimated_states(:,k));
        second_derivative_of_f4 = d2_dx2(f4,pseudotrue_estimated_states(:,k));

        first_derivative_of_h = d_dx(h,pseudotrue_estimated_states(:,k));
        second_derivative_of_h1 = d2_dx2(h1,pseudotrue_estimated_states(:,k));
        second_derivative_of_h2 = d2_dx2(h2,pseudotrue_estimated_states(:,k));
    else
        % use the analytical derivative expressions
        first_derivative_of_f = df_dx(pseudotrue_estimated_states(:,k));

        second_derivative_of_f1 = df2_dx2(pseudotrue_estimated_states(:,k));
        second_derivative_of_f2 = df2_dx2(pseudotrue_estimated_states(:,k));
        second_derivative_of_f3 = df2_dx2(pseudotrue_estimated_states(:,k));
        second_derivative_of_f4 = df2_dx2(pseudotrue_estimated_states(:,k));

        first_derivative_of_h = dh_dx(pseudotrue_estimated_states(:,k));

        second_derivative_of_h1 = d2h1_dx2(pseudotrue_estimated_states(:,k));
        second_derivative_of_h2 = d2h2_dx2(pseudotrue_estimated_states(:,k));
    end

    % calculate k+1 terms
    if k < simulation_length
        first_elementf = (pseudotrue_estimated_states(:,k+1) - f(pseudotrue_estimated_states(:,k)))'*invQ*[1;0;0;0]*second_derivative_of_f1;
        second_elementf = (pseudotrue_estimated_states(:,k+1) - f(pseudotrue_estimated_states(:,k)))'*invQ*[0;1;0;0]*second_derivative_of_f2;
        third_elementf = (pseudotrue_estimated_states(:,k+1) - f(pseudotrue_estimated_states(:,k)))'*invQ*[0;0;1;0]*second_derivative_of_f3;
        fourth_elementf = (pseudotrue_estimated_states(:,k+1) - f(pseudotrue_estimated_states(:,k)))'*invQ*[0;0;0;1]*second_derivative_of_f4;
        second_derivative_sum_of_f = first_elementf + second_elementf + third_elementf + fourth_elementf;
    end

    first_elementh = (true_measurements_in_polar(:,k) - h(pseudotrue_estimated_states(:,k)))'*invR*[1;0]*second_derivative_of_h1;
    second_elementh = (true_measurements_in_polar(:,k) - h(pseudotrue_estimated_states(:,k)))'*invR*[0;1]*second_derivative_of_h2;
    second_derivative_sum_of_h = first_elementh + second_elementh;

    % no k-1 term
    if k == 1
        Amatrixkthblock = (- first_derivative_of_f' *invQ * first_derivative_of_f)...
            +second_derivative_sum_of_f...
            -first_derivative_of_h' *invR *first_derivative_of_h...
            +second_derivative_sum_of_h...
            -inv(P0);

        % no k+1 terms
    elseif k ==simulation_length
        Amatrixkthblock = -invQ ...
            -first_derivative_of_h' *invR *first_derivative_of_h...
            +second_derivative_sum_of_h;

        % general case
    else
        Amatrixkthblock = -invQ + (- first_derivative_of_f' *invQ * first_derivative_of_f)...
            +second_derivative_sum_of_f...
            -first_derivative_of_h' *invR *first_derivative_of_h...
            +second_derivative_sum_of_h;

    end

    row_start = (k-1) * state_dimension + 1;
    row_end = k * state_dimension;
    col_start = row_start;
    col_end = row_end;

    Amatrix(row_start:row_end, col_start:col_end) = Amatrixkthblock;

    % sub and super diagonal terms calculation
    if k < simulation_length
        Amatrixkthblocknondiagonalterm = first_derivative_of_f'*invQ;
        super_row_start = row_start;
        super_row_end = row_end;
        super_col_start = col_end + 1;
        super_col_end = col_end + state_dimension;
        Amatrix(super_row_start:super_row_end, super_col_start:super_col_end) = Amatrixkthblocknondiagonalterm;
        Amatrix(super_col_start:super_col_end,super_row_start:super_row_end) = Amatrixkthblocknondiagonalterm';
    end


    if k == 1
        b_element = ...
            +first_derivative_of_f'*invQ*(pseudotrue_estimated_states(:,k+1)-f(pseudotrue_estimated_states(:,k)))...
            +first_derivative_of_h'*invR*(true_measurements_in_polar(:,k)-h(pseudotrue_estimated_states(:,k)))...
            -inv(P0)*(pseudotrue_estimated_states(:,k) - x0);

    elseif k == simulation_length
        b_element = -invQ * (pseudotrue_estimated_states(:,k)-f(pseudotrue_estimated_states(:,k-1)))...
            +first_derivative_of_h'*invR*(true_measurements_in_polar(:,k)-h(pseudotrue_estimated_states(:,k)));

    else
        b_element = -invQ * (pseudotrue_estimated_states(:,k)-f(pseudotrue_estimated_states(:,k-1)))...
            +first_derivative_of_f'*invQ*(pseudotrue_estimated_states(:,k+1)-f(pseudotrue_estimated_states(:,k)))...
            +first_derivative_of_h'*invR*(true_measurements_in_polar(:,k)-h(pseudotrue_estimated_states(:,k)));
    end
    bvector(row_start:row_end) = b_element;
    Bmatrix(row_start:row_end, col_start:col_end) = first_derivative_of_h'*invR*Rbar*invR*first_derivative_of_h;

end
Bmatrix=Bmatrix+bvector*bvector';

if apply_tridiagonal_method
    % calculate the inverse of A using the symmetric block tridiagonal property of A (Meurant, 1992)
    Ainv = calculate_Ainv(Amatrix, state_dimension);
else
    % calculate the inverse of A using the built-in inverse function
    Ainv = inv(Amatrix);
end

MCRLB = Ainv * Bmatrix * Ainv;

end