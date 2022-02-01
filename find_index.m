function [total_index] = find_index(L,m,EqS)

    LL = L + m -1;
    if (mod(LL-m,2) == EqS)
        total_index =( 5*(LL-m+2-EqS)/2 );
    else
        total_index =( 5*(LL-m+1-EqS)/2 );
    end

end