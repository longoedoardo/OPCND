function x = mono_next_grlex(m, x)
    % Validate inputs fast
    if m < 1
        error('MONO_NEXT_GRLEX: M must be >= 1.');
    end
    if any(x < 0)
        error('MONO_NEXT_GRLEX: All components of x must be >= 0.');
    end

    % Find index of rightmost nonzero entry (built-in find is efficient)
    i = find(x > 0, 1, 'last');

    if isempty(i)
        % x is zero vector -> next is unit in last position
        x(m) = 1;
        return
    end

    if i == 1
        t = x(1) + 1;
        im1 = m;
    else
        t = x(i);
        im1 = i - 1;
    end

    x(i) = 0;
    x(im1) = x(im1) + 1;
    x(m) = x(m) + t - 1;
end
